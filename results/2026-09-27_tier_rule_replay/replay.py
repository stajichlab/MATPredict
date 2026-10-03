"""Replay the proposed allele-absent tier rule on existing detection reports.

Run once per (run, code tree) with PYTHONPATH=<tree>/src so the family roster
and expected_genes_for_idiomorph come from the SAME tree that produced the run.

Reimplements tiering.assign_tier plus the post-tier caps in pipeline._build /
flank_carried from report fields only:
  base tier (core_found / isolated single hit / any_gene_unpolished / flank)
  -> fragmented downgrade
  -> idiomorph-margin cap (no classifier verdict, margin < family.min_idiomorph_margin)
  -> polished_genes == 0 -> low
  -> partial_locus (strength) -> cap medium
  -> relaxed pass -> cap medium
  -> idiomorph_unmodelled (flank-carried) -> low

"unpolished" = gene_evidence status == "unpolished" (polish.STATUS_UNPOLISHED).

Variants (applied only when the call's idiomorph is a single named idiomorph):
  A  the rule as proposed: ignore UNPOLISHED genes whose present_in_idiomorphs
     excludes the called idiomorph when deciding any_gene_unpolished.
  B  A, plus drop allele-absent genes (any status) from the found set used to
     compute the expected core (expected_genes_for_idiomorph).
  B' A, plus drop only WEAK (<50% identity) allele-absent genes from that set.

Usage: replay.py RUN_LABEL DB_ROOT REPORT_GLOB_OR_TARBALL OUT_TSV
"""
import csv, glob, io, sys, tarfile
import yaml
from pathlib import Path

from MATPredict.detect.family_registry import load_all_families, expected_genes_for_idiomorph

label, db_root, src, out = sys.argv[1:5]
fams = {}
for f in load_all_families(Path(db_root)):
    fams[f"{f.key.phylum}:{f.key.locus_name}" if hasattr(f.key, "phylum") else str(f.key)] = f
TIER_DOWN = {"high": "medium", "medium": "low", "low": "low"}


def fam_for(key):
    if key in fams:
        return fams[key]
    for k, f in fams.items():
        if k.endswith(":" + key.split(":")[-1]) and k.split(":")[0] == key.split(":")[0]:
            return f
    return None


def reports():
    if src.endswith(".tar.zst"):
        import subprocess
        data = subprocess.run(["zstd", "-dc", src], capture_output=True, check=True).stdout
        with tarfile.open(fileobj=io.BytesIO(data)) as tf:
            for m in tf.getmembers():
                if m.name.endswith("detection_report.yaml"):
                    yield m.name, yaml.safe_load(tf.extractfile(m))
    else:
        for p in sorted(glob.glob(src, recursive=True)):
            yield p, yaml.safe_load(open(p))


def base_tier(fam, found, not_searchable, unpolished_any):
    exp = expected_genes_for_idiomorph(fam, found)
    expected_core = {g["name"] for g in exp if g["role"] == "core_MAT" and not g.get("optional")}
    core_genes = expected_core - set(not_searchable)
    relaxed_core = core_genes != expected_core
    core_found = core_genes.issubset(set(found))
    isolated = len(found) <= 1 and len(fam.genes) > 1
    if not core_found:
        # fraction_found == 0 never reaches a report; isolated -> low
        return ("low" if isolated else "medium"), "core_not_found"
    if isolated and relaxed_core:
        return "low", "isolated"
    if unpolished_any:
        return "medium", "unpolished"
    unsearch = set(not_searchable)
    flank_req = any(g["role"] == "flanking_conserved" and g["name"] not in unsearch
                    and not g.get("optional") for g in fam.genes)
    if flank_req:
        ok = any(g["role"] == "flanking_conserved" and g["name"] in found for g in fam.genes)
        return ("high", "flank_ok") if ok else ("medium", "flank_missing")
    return "high", "no_flank_req"


def final_tier(fam, x, found, unpolished_any):
    t, why = base_tier(fam, found, x.get("genes_not_searchable") or [], unpolished_any)
    if x.get("fragmented") and t != "low":
        t = TIER_DOWN[t]; why += "+fragmented"
    m = x.get("idiomorph_margin")
    if not x.get("idiomorph_classifier") and m is not None and t == "high":
        if m < getattr(fam, "min_idiomorph_margin", 0) or 0:
            t = "medium"; why += "+margin"
    if x.get("polished_genes") == 0:
        t = "low"; why += "+zero_modelled"
    if x.get("locus_class") == "partial_locus" and t == "high":
        t = "medium"; why += "+partial"
    if x.get("detection_pass") == "relaxed" and t == "high":
        t = "medium"; why += "+relaxed"
    if x.get("idiomorph_unmodelled"):
        t = "low"; why += "+flank_carried"
    return t, why


rows = []
for path, r in reports():
    genome = path.split("/")[-2]
    for x in r.get("detected") or []:
        fam = fam_for(x["family"])
        if fam is None:
            rows.append(dict(run=label, genome=genome, family=x["family"], note="family_not_in_roster"))
            continue
        found = list(x.get("genes_found") or [])
        ev = x.get("gene_evidence") or []
        unpol = {e["gene"] for e in ev if e.get("status") == "unpolished"}
        idio = x.get("idiomorph")
        pin = {g["name"]: set(g.get("present_in_idiomorphs") or []) for g in fam.genes}
        absent = set()
        if idio and idio not in ("undetermined",) and "+" not in str(idio):
            absent = {g for g in set(found) | {e["gene"] for e in ev}
                      if pin.get(g) and idio not in pin[g]}
        # model_idiomorph_alternatives families already exclude the losing
        # (other-allele) half from the unpolished check (pipeline live_genes).
        unpol_cur = unpol - absent if getattr(fam, "model_idiomorph_alternatives", False) else unpol
        cur, why = final_tier(fam, x, found, bool(unpol_cur))
        a, wa = final_tier(fam, x, found, bool(unpol - absent))
        b, wb = final_tier(fam, x, [g for g in found if g not in absent], bool(unpol - absent))
        weak_absent = {e["gene"] for e in ev if e["gene"] in absent
                       and e.get("identity") is not None and float(e["identity"]) < 50}
        bp, wbp = final_tier(fam, x, [g for g in found if g not in weak_absent], bool(unpol - absent))
        ign = [(e["gene"], e.get("identity"), e.get("status")) for e in ev if e["gene"] in absent]
        rows.append(dict(
            run=label, genome=genome, family=x["family"], idiomorph=idio,
            reported=x["confidence"], reproduced=cur, reproduced_ok=cur == x["confidence"],
            why_current=why, variant_A=a, why_A=wa, variant_B=b, why_B=wb, variant_Bp=bp, why_Bp=wbp,
            absent_genes=";".join(f"{g}:{i}:{s}" for g, i, s in ign),
            max_absent_identity=max([i for _, i, _ in ign if i is not None], default=""),
            unpolished=";".join(sorted(unpol)), polished_genes=x.get("polished_genes"),
            locus_class=x.get("locus_class"), contig=x.get("contig"), start=x.get("start"),
            end=x.get("end"),
            core_evidence=";".join(f"{e['gene']}:{e.get('identity')}:{e.get('status')}"
                                   for e in ev if e.get("role") == "core_MAT"),
            flanks=";".join(sorted(e["gene"] for e in ev if str(e.get("role", "")).startswith("flanking"))),
        ))
keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("run", "genome", "family"), k))
with open(out, "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=keys, delimiter="\t"); w.writeheader(); w.writerows(rows)
print(label, len(rows), "calls")
