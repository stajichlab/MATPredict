"""Replay the IMPLEMENTED allele-absent tier rule on existing reports.

Unlike results/2026-09-27_tier_rule_replay/replay.py, the base tier and the
ignore set come from the real functions in the new code
(`tiering.assign_tier`, `tiering.allele_absent_genes_to_ignore`), run with
PYTHONPATH=<polish-scope-cuts>/src. The family roster is the one from the
tree that produced each run (DB_ROOT). The post-tier caps of
`pipeline._build` / `flank_carried` (unchanged by this rule) are re-applied
from report fields, as in the earlier replay.

"modelled" = gene_evidence status polished_agree/_disagree/_single, or method
diamond_proteome (pipeline._modelled_gene_names).

Usage: replay_real.py RUN_LABEL DB_ROOT REPORT_GLOB_OR_TARBALL OUT_TSV
"""
import csv, glob, io, subprocess, sys, tarfile
from pathlib import Path

import yaml

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import load_all_families
from MATPredict.detect.scoring import FamilyScore
from MATPredict.detect.tiering import allele_absent_genes_to_ignore, assign_tier

label, db_root, src, out = sys.argv[1:5]
fams = {f"{f.key.phylum}:{f.key.locus_name}": f for f in load_all_families(Path(db_root))}
TIER_DOWN = {"high": "medium", "medium": "low", "low": "low"}
MODELLED = {"polished_agree", "polished_disagree", "polished_single"}


def reports():
    if src.endswith(".tar.zst"):
        data = subprocess.run(["zstd", "-dc", src], capture_output=True, check=True).stdout
        with tarfile.open(fileobj=io.BytesIO(data)) as tf:
            for m in tf.getmembers():
                if m.name.endswith("detection_report.yaml"):
                    yield m.name, yaml.safe_load(tf.extractfile(m))
    else:
        for p in sorted(glob.glob(src, recursive=True)):
            yield p, yaml.safe_load(open(p))


def caps(fam, x, t):
    """Post-tier caps from pipeline._build / flank_carried, unchanged."""
    m = x.get("idiomorph_margin")
    if not x.get("idiomorph_classifier") and m is not None and t == "high":
        if m < (getattr(fam, "min_idiomorph_margin", 0) or 0):
            t = "medium"
    if x.get("polished_genes") == 0:
        t = "low"
    if x.get("locus_class") == "partial_locus" and t == "high":
        t = "medium"
    if x.get("detection_pass") == "relaxed" and t == "high":
        t = "medium"
    if x.get("idiomorph_unmodelled"):
        t = "low"
    return t


def tier(fam, x, unpolished, ignore):
    score = FamilyScore(fam.key, float(x.get("fraction_found") or 1.0),
                        list(x.get("genes_found") or []), list(x.get("genes_missing") or []),
                        list(x.get("genes_not_searchable") or []))
    t = assign_tier(score, fam, GeneCluster(x.get("contig", "c"), 1, 2, []),
                    any_gene_unpolished=bool(unpolished - ignore),
                    fragmented=bool(x.get("fragmented")), ignore_genes=ignore)
    return caps(fam, x, t)


rows = []
for path, r in reports():
    genome = path.split("/")[-2]
    for x in r.get("detected") or []:
        fam = fams.get(x["family"])
        if fam is None:
            continue
        ev = x.get("gene_evidence") or []
        idio = x.get("idiomorph")
        pin = {g["name"]: set(g.get("present_in_idiomorphs") or []) for g in fam.genes}
        unpol = {e["gene"] for e in ev if e.get("status") == "unpolished"}
        if getattr(fam, "model_idiomorph_alternatives", False) and idio:
            # pipeline live_genes: the losing half of the pair is not counted
            unpol = {g for g in unpol if not (pin.get(g) and idio not in pin[g])}
        ident = {}
        for e in ev:
            i = e.get("identity")
            if i is not None:
                ident[e["gene"]] = max(float(i), ident.get(e["gene"]) or 0.0)
            else:
                ident.setdefault(e["gene"], None)
        modelled = {e["gene"] for e in ev
                    if e.get("status") in MODELLED or e.get("method") == "diamond_proteome"}
        ignore = allele_absent_genes_to_ignore(fam, idio, ident, modelled, locus_class=x.get("locus_class"))
        cur = tier(fam, x, unpol, frozenset())
        new = tier(fam, x, unpol, ignore)
        rows.append(dict(
            run=label, genome=genome, family=x["family"], idiomorph=idio,
            reported=x["confidence"], reproduced=cur, reproduced_ok=cur == x["confidence"],
            new=new, ignored=";".join(f"{g}:{ident.get(g)}" for g in sorted(ignore)),
            locus_class=x.get("locus_class"), polished_genes=x.get("polished_genes"),
            contig=x.get("contig"), start=x.get("start"), end=x.get("end"),
            core_evidence=";".join(f"{e['gene']}:{e.get('identity')}:{e.get('status')}"
                                   for e in ev if e.get("role") == "core_MAT"),
            flanks=";".join(sorted(f"{e['gene']}:{e.get('identity')}:{e.get('status')}"
                                   for e in ev if str(e.get("role", "")).startswith("flanking"))),
        ))
keys = list(rows[0]) if rows else ["run"]
with open(out, "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=keys, delimiter="\t"); w.writeheader(); w.writerows(rows)
ok = sum(r["reproduced_ok"] for r in rows)
print(f"{label}: {len(rows)} calls, reproduced {ok}/{len(rows)}, "
      f"risers {sum(r['reproduced_ok'] and r['reported'] == 'medium' and r['new'] == 'high' for r in rows)}")
