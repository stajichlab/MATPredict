"""Fraction of each shipped profile (alpha box, HMG box) covered by the best alignment of a protein: for the curated MAT1-1-1 / MAT1-2-1
records and for the alpha/HMG models of the genomes re-classed in the first patched run.  python cov_check.py DB REPORTS_DIR LIB OUT
(run on the HPCC, needs pyhmmer)."""
import gzip, shutil, sys, tempfile
from pathlib import Path
import pyhmmer, yaml
from MATPredict.detect import partial_core, classifier_build as cb
from MATPredict.detect.family_registry import FamilyKey, load_all_families
from MATPredict.detect.pipeline import _translate_model
from MATPredict.detect.reference_fasta import build_reference_fasta

db, runs, lib, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
family = next(f for f in load_all_families(db) if f.key == FamilyKey("Ascomycota", "MAT"))
clf = partial_core.load_core_profiles(db, family)
hm = {gene: h for _i, (gene, h) in clf.models.items()}


def cover(hmm, protein):
    al = pyhmmer.easel.Alphabet.amino()
    blk = pyhmmer.easel.DigitalSequenceBlock(al, [pyhmmer.easel.TextSequence(name=b"p", sequence=protein).digitize(al)])
    pl = pyhmmer.plan7.Pipeline(al, Z=1, E=1e9, domE=1e9, bias_filter=False, F1=1.0, F2=1.0, F3=1.0)
    best = (0.0, 0.0)
    for hit in pl.search_hmm(hmm, blk):
        for d in hit.domains:
            c = (d.alignment.hmm_to - d.alignment.hmm_from + 1) / hmm.M
            if d.score > best[0]:
                best = (float(d.score), c)
    return best


rows = ["kind\tid\tgene\tprotein_aa\talpha_bits\talpha_cov\thmg_bits\thmg_cov"]
ref = Path(tempfile.mkdtemp()) / "r.faa"
build_reference_fasta(db, ref, family_keys={family.key})
for h, s in cb.read_fasta(ref).items():
    g = h.split("|")[-1]
    if g in ("MAT1-1-1", "MAT1-2-1"):
        a, hg = cover(hm["MAT1-1-1"], s), cover(hm["MAT1-2-1"], s)
        rows.append(f"curated\t{h.split('|')[0]}\t{g}\t{len(s)}\t{a[0]:.1f}\t{a[1]:.2f}\t{hg[0]:.1f}\t{hg[1]:.2f}")
tmp = Path(tempfile.mkdtemp())
for rep in sorted(runs.glob("*/detection_report.yaml")):
    g = rep.parent.name
    d = yaml.safe_load(open(rep)) if False else yaml.load(open(rep), Loader=getattr(yaml, "CSafeLoader", yaml.SafeLoader))
    if not any(x.get("homothallic_note") for x in d.get("detected") or []):
        continue
    fa = tmp / "g.fna"
    with gzip.open(lib / f"{g}.fa.gz", "rt") as fi, open(fa, "w") as fo:
        shutil.copyfileobj(fi, fo)
    cache = {}
    for x in d["detected"]:
        n = x.get("homothallic_note")
        if not n:
            continue
        for key in ("alpha_box_gene", "hmg_box_gene"):
            name = n[key]["gene"]
            e = next(e for e in x["gene_evidence"] if e["gene"] == name and e.get("exons"))
            ex = [(a["start"], a["end"]) for a in e["exons"]]
            span = type("S", (), dict(contig=e["contig"], start=e["start"], end=e["end"], strand=e["strand"], frameshifts=0,
                                      cds_blocks=(), exons=[type("E", (), dict(start=a, end=b)) for a, b in ex]))
            prot = _translate_model(fa, span, d["genetic_code"] or 1, cache)
            a, hg = cover(hm["MAT1-1-1"], prot), cover(hm["MAT1-2-1"], prot)
            rows.append(f"patched\t{g}\t{name}\t{len(prot)}\t{a[0]:.1f}\t{a[1]:.2f}\t{hg[0]:.1f}\t{hg[1]:.2f}")
open(out, "w").write("\n".join(rows) + "\n")
print(len(rows))
