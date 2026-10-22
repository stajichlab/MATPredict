"""Per-model metrics behind `homothallic.evaluate` for every locus of the patched run: which test each core model passes or fails.
Run on the HPCC (needs pyhmmer): python why_not.py RUNS_DIR DB_ROOT LIB OUT_TSV [genome list file]
One row per modelled core gene of every detected Ascomycota:MAT locus: protein length, coding length, ORF fraction, bits on the
alpha and HMG profiles, and the first failing condition (`ok`, `polish_status`, `short`, `orf_truncated`, `profile_floor`,
`profile_margin`)."""
import csv, gzip, sys
from pathlib import Path
import yaml

from MATPredict.detect import homothallic as hm, partial_core
from MATPredict.detect.classifier import score_proteins
from MATPredict.detect.family_registry import FamilyKey, load_all_families
from MATPredict.detect.pipeline import _translate_model

runs, db, lib, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
only = set(open(sys.argv[5]).read().split()) if len(sys.argv) > 5 else None
L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
family = next(f for f in load_all_families(db) if f.key == FamilyKey("Ascomycota", "MAT"))
clf = partial_core.load_core_profiles(db, family)
cls_of = {i: next(g.get("gene_class") for g in family.genes if g["name"] == gene) for i, (gene, _h) in clf.models.items()}
import tempfile, shutil, os
rows = []
tmp = Path(tempfile.mkdtemp())
for rep in sorted(runs.glob("*/detection_report.yaml")):
    g = rep.parent.name
    if only is not None and g not in only:
        continue
    d = yaml.load(open(rep), Loader=L)
    locs = [x for x in d.get("detected") or [] if x["family"] == "Ascomycota:MAT"]
    if not locs:
        continue
    fa = tmp / "g.fna"
    with gzip.open(lib / f"{g}.fa.gz", "rt") as fin, open(fa, "w") as fout:
        shutil.copyfileobj(fin, fout)
    cache = {}
    for x in locs:
        for e in x.get("gene_evidence") or []:
            if e["role"] != "core_MAT":
                continue
            row = dict(genome=g, contig=x["contig"], start=x["start"], gene=e["gene"], status=e["status"], lost_to=e.get("lost_to"),
                       locus_class=x["locus_class"], reclassed=bool(x.get("homothallic_note")))
            if e["status"] not in hm.MODELLED_STATUSES or not e.get("exons"):
                row.update(fail="polish_status"); rows.append(row); continue
            ex = [(a["start"], a["end"]) for a in e["exons"]]
            span = type("S", (), dict(contig=e["contig"], start=e["start"], end=e["end"], strand=e["strand"], frameshifts=0,
                                      cds_blocks=(), exons=[type("E", (), dict(start=a, end=b)) for a, b in ex]))
            prot = _translate_model(fa, span, d["genetic_code"] or 1, cache)
            sc = score_proteins(clf, [prot]) if prot else {}
            s = {cls_of[i]: v for i, v in sc.items() if cls_of.get(i)}
            cds = sum(b - a + 1 for a, b in ex) // 3
            m = hm.CoreModel(e["gene"], e["contig"], e["start"], e["end"], len(prot), cds, s)
            if len(prot) < hm.MIN_PROTEIN_AA:
                fail = "short"
            elif len(prot) < hm.MIN_ORF_FRACTION * cds:
                fail = "orf_truncated"
            else:
                c = hm.confirmed_class(m)
                if c:
                    fail = "ok"
                else:
                    best = max(s, key=s.get) if s else None
                    fail = "profile_floor" if best is None or s[best] < hm.MIN_PROFILE_BITS[best] else "profile_margin"
            row.update(protein_aa=len(prot), cds_aa=cds, alpha=round(s.get("alpha_box", 0.0), 1), hmg=round(s.get("HMG_box", 0.0), 1),
                       strand=e["strand"], fail=fail, confirmed=hm.confirmed_class(m))
            rows.append(row)
cols = []
for r in rows:
    for k in r:
        if k not in cols:
            cols.append(k)
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t", restval="")
    w.writeheader(); w.writerows(rows)
print(len(rows), "core models")
