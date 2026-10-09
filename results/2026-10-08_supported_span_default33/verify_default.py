#!/usr/bin/env python3
"""The new default (33 bits, no flag) against (a) the earlier explicit --supported-min-bitscore 33 runs, field by field, and
(b) current main (calls, cluster spans, genes, merged_from, withheld loci), plus containment on the default runs.
usage: verify_default.py (in this folder: basidio/ other/; earlier runs in ../2026-10-08_supported_span/floor33 and
../2026-10-08_supported_span_other_groups/{floor33,base}/ )"""
import glob, yaml
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/"
SETS = {
    "Basidiomycota panel (34)": ("basidio", R + "2026-10-08_supported_span/floor33", R + "2026-10-08_supported_span_withheld/newdb"),
    "Ascomycota + Mucoromycota (110)": ("other", R + "2026-10-08_supported_span_other_groups/floor33", R + "2026-10-08_supported_span_other_groups/base"),
}
def load(d): return {f.split("/")[-2]: yaml.safe_load(open(f)) for f in glob.glob(f"{d}/runs/*/detection_report.yaml")}
def strip(r):   # drop run provenance (timing, host) before comparing
    r = dict(r); r.pop("run", None); return r
def callkey(x): return (x["family"], x["contig"], x["start"], x["end"], x["confidence"], x["idiomorph"], x.get("locus_class"),
                        tuple(sorted(x.get("genes_found") or [])), tuple(sorted(str(m) for m in (x.get("merged_from") or []))))
def sig(r): return ([callkey(x) for x in r.get("detected") or []],
                    sorted((s["family"], s["contig"], s["start"], s["end"], s.get("withheld_reason"), tuple(s.get("genes_found") or [])) for s in r.get("suppressed_loci") or []))
for name, (d, f33, base) in SETS.items():
    D, E, B = load(d), load(f33), load(base)
    print(f"\n{name}: default-run reports {len(D)}, explicit-33 reports {len(E)}, main/earlier-baseline {len(B)}")
    same = [g for g in D if g in E and strip(D[g]) == strip(E[g])]
    print(f"   default vs explicit floor 33: {len(same)} of {len(D)} reports identical in every field (supported_span included)")
    common = [g for g in D if g in B]
    bad = [g for g in common if sig(D[g]) != sig(B[g])]
    print(f"   calls, cluster spans, genes, merged_from and withheld loci vs {base.split('/')[-2]}/{base.split('/')[-1]}: {len(common)} genomes compared, {len(bad)} differ {bad[:3]}")
    n = out = 0; floors = set()
    for g, r in D.items():
        for x in (r.get("detected") or []) + (r.get("suppressed_loci") or []):
            sp = x.get("supported_span")
            if sp: floors.add(sp["min_bitscore"])
            for e in x.get("gene_evidence") or []:
                if sp and e.get("contig") == x["contig"]:
                    n += 1; out += not (sp["start"] <= e["start"] and e["end"] <= sp["end"])
    print(f"   containment: {n} reported genes checked, {out} outside the supported span; floors recorded in the reports: {sorted(floors)}")
