#!/usr/bin/env python3
"""Read-based MTL zygosity for C. albicans isolates vs the assembly call.

Depth per region comes from map2.slurm (depth2.tsv: 3M randomly sampled R1
reads per run, minimap2 -ax sr, MAPQ>=20, primary only). depth.tsv (map.slurm,
first 2M reads) is used only for isolates absent from depth2.tsv and only when
its control depth is not inflated (see NOTE.md).

Each idiomorph is scored as the mean of its two blocks (MTL genes + idiomorph
PAP/OBP/PIK alleles), divided by the chr1 control depth.
  present : ratio >= 0.20
  absent  : ratio <  0.05
  between : unresolved
"""
import csv, collections, sys

A = ("a_MTLa2_MTLa1", "a_PAPa_OBPa_PIKa")
AL = ("alpha_MTLalpha2", "alpha_MTLalpha1", "alpha_OBP_PIK_PAP")

def load(path):
    d = collections.defaultdict(dict); meta = {}
    for acc, geno, run, region, depth in csv.reader(open(path), delimiter="\t"):
        d[acc][region] = float(depth); meta[acc] = (geno, run)
    return d, meta

d2, m2 = load("depth2.tsv")
d1, m1 = load("depth.tsv")
rows = []
for acc in sorted(set(d2) | set(d1)):
    if acc in d2:
        dep, (geno, run), src = d2[acc], m2[acc], "round2"
    else:
        dep, (geno, run), src = d1[acc], m1[acc], "round1"
        if dep["control"] > 40:          # position-ordered submission; unusable
            continue
    c = dep["control"]
    ra = sum(dep[r] for r in A) / len(A) / c
    ral = sum(dep[r] for r in AL) / len(AL) / c
    call = lambda r: "present" if r >= 0.20 else ("absent" if r < 0.05 else "unresolved")
    ca, cal = call(ra), call(ral)
    reads = {("present", "present"): "a/alpha", ("present", "absent"): "a/a",
             ("absent", "present"): "alpha/alpha"}.get((ca, cal), "unresolved")
    rows.append((acc, run, src, geno, round(c, 2), round(ra, 2), round(ral, 2), reads))

w = csv.writer(sys.stdout, delimiter="\t")
w.writerow(["assembly", "run", "depth_source", "assembly_call", "control_depth",
            "a_ratio", "alpha_ratio", "read_zygosity"])
w.writerows(rows)
conf = collections.Counter((r[3], r[7]) for r in rows)
print("\n# confusion: assembly_call x read_zygosity", file=sys.stderr)
for (g, z), n in sorted(conf.items()):
    print(f"#  {g:8s} {z:12s} {n}", file=sys.stderr)
