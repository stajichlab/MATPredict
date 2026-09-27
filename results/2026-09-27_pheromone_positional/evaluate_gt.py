#!/usr/bin/env python3
"""Label ground-truth STE3 loci (mating vs non-mating) and compare precursor
counts near them with chance (random windows of equal size)."""
import csv
import glob
import os
import sys
from collections import defaultdict

D = sys.argv[1] if len(sys.argv) > 1 else "out_gt"
# mating receptors: record coordinates (same assembly) or >=95% to the strain's own record receptor
COORD = {
    "GCA_016772295.1_ASM1677229v1": [("JAAGWA010000010.1", 1806154, 1823946)],
    "GCF_000143185.2_Schco3": [("NW_026089548.1", 218046, 220353), ("NW_026089548.1", 246348, 248430)],
}
OWNREF = {
    "GCF_000091045.1_ASM9104v1": "REF|40410_jec21_MAT_alpha|",
    "GCF_000328475.2_Umaydis521_2.0": "REF|5270_521_aLocus_a1|",
    "GCA_921037615.3_Hybrid_genome_assembly_and_annotation_02": "REF|5286_cbs-14_redPR_A1|",
    "GCA_000988875.2_ASM98887v2": "REF|5286_nbrc-0880_redPR_A2|",
}
CLS = ["T", "T2", "R", "T2|R", "L", "Hx"]
WS = ("10kb", "20kb", "50kb")


def pos(r, c, W):
    if c == "R":
        return int(r[f"R_{W}"]) >= 2
    if c == "T2|R":
        return int(r[f"T2_{W}"]) > 0 or int(r[f"R_{W}"]) >= 2
    return int(r[f"{c}_{W}"]) > 0


def label(r):
    a = r["asm"]
    if a in COORD:
        s, e = int(r["start"]), int(r["end"])
        return any(r["contig"] == c and s <= ce and e >= cs for c, cs, ce in COORD[a])
    return r["best_ref"].startswith(OWNREF[a]) and float(r["best_ref_ident"]) >= 0.95


rows_out = []
summary = defaultdict(lambda: defaultdict(int))
for f in sorted(glob.glob(os.path.join(D, "*.loci.tsv"))):
    asm = os.path.basename(f).replace(".loci.tsv", "")
    loci = list(csv.DictReader(open(f), delimiter="\t"))
    rnd = list(csv.DictReader(open(f.replace(".loci.tsv", ".random.tsv")), delimiter="\t"))
    for r in loci:
        m = label(r)
        r["mating"] = "mating" if m else "other"
        rows_out.append(r)
        for W in WS:
            for c in CLS:
                summary[(asm, r["mating"], W, c)]["n"] += 1
                summary[(asm, r["mating"], W, c)]["pos"] += pos(r, c, W)
    for W in WS:
        for c in CLS:
            summary[(asm, "random", W, c)]["n"] += len(rnd)
            summary[(asm, "random", W, c)]["pos"] += sum(pos(x, c, W) for x in rnd)

with open(os.path.join(D, "..", "gt_loci_labelled.tsv"), "w") as fo:
    keys = ["asm", "contig", "strand", "start", "end", "mating", "best_query", "best_ident", "best_ref",
            "best_ref_ident"] + [f"{c}_{w}" for w in WS for c in ["T", "T2", "L", "R", "H1", "H01", "Hx"]]
    fo.write("\t".join(keys) + "\n")
    for r in rows_out:
        fo.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")

print("asm\tgroup\tW\t" + "\t".join(CLS))
tot = defaultdict(lambda: defaultdict(lambda: [0, 0]))
for asm in sorted({k[0] for k in summary}):
    for g in ("mating", "other", "random"):
        for W in WS:
            cells = []
            for c in CLS:
                d = summary.get((asm, g, W, c))
                if not d or not d["n"]:
                    cells.append("-")
                    continue
                cells.append(f"{d['pos']}/{d['n']}")
                tot[(g, W)][c][0] += d["pos"]
                tot[(g, W)][c][1] += d["n"]
            print(f"{asm[:28]}\t{g}\t{W}\t" + "\t".join(cells))
print("\nPOOLED (random = per-window rate)")
for g in ("mating", "other", "random"):
    for W in WS:
        print(f"{g}\t{W}\t" + "\t".join(f"{tot[(g,W)][c][0]}/{tot[(g,W)][c][1]} ({100*tot[(g,W)][c][0]/max(1,tot[(g,W)][c][1]):.1f}%)" for c in CLS))
