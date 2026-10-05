#!/usr/bin/env python3
"""Labelled-panel test of the strict-CAAX rule (T, +-10 kb) on STE3-like loci.

Inputs:
  ../2026-09-27_pheromone_positional/gt_loci_labelled.tsv   old panel, 6 genomes, 34 loci
  out_new/<asm>.loci.tsv                                     4 new Agaricomycete genomes
Label (new genomes): mating = locus overlaps the curated record's receptor cluster.
Everything else in the genome is "other" (assumed non-mating; see README).

Independence of positives: the new records (Trametes, Grifola, Russula,
Heterobasidion) were curated from CAAX-positional evidence, so their mating
receptors are NOT independent of the rule. They are reported but kept out of
the sensitivity estimate. Their "other" copies are used as negatives.
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(sys.argv[0]))
OLD = os.path.join(HERE, "..", "2026-09-27_pheromone_positional", "gt_loci_labelled.tsv")
NEW_DIR = os.path.join(HERE, "out_new")

# record receptor-cluster coordinates (db/Basidiomycota/*/<record>/metadata.yaml)
NEW_COORD = {
    # GCF assemblies use RefSeq contig names; the record sequence was matched to the
    # genome (record base 1 = 1,556,459 / 825,882), so only the contig name changes.
    "GCF_000271585.1_Trametes_versicolor_v1.0": [("NW_007360328.1", 1556459, 1581006)],
    "GCA_001683735.1_ASM168373v1": [("LUGG01000005.1", 512153, 547096)],
    "GCA_984573805.1_gfRusNobi1.hap1.1": [("OZ475200.1", 1514818, 1529433)],
    "GCF_000320585.1_Heterobasidion_irregulare_v2.0": [("NW_009258203.1", 825882, 849398)],
}
# old panel: which genomes are Agaricomycete (the family the scan is enabled for)
OLD_AGARICO = {"GCA_016772295.1_ASM1677229v1", "GCF_000143185.2_Schco3"}


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def rule(r, w="10kb"):
    return int(r[f"T_{w}"]) > 0


rows = []
for r in csv.DictReader(open(OLD), delimiter="\t"):
    r["set"] = "old"
    r["agarico"] = r["asm"] in OLD_AGARICO
    r["independent"] = True
    rows.append(r)
for asm, coords in NEW_COORD.items():
    f = os.path.join(NEW_DIR, asm + ".loci.tsv")
    if not os.path.exists(f):
        print(f"missing {f}", file=sys.stderr)
        continue
    for r in csv.DictReader(open(f), delimiter="\t"):
        r["asm"] = asm
        s, e = int(r["start"]), int(r["end"])
        r["mating"] = "mating" if any(
            r["contig"] == c and s <= ce and e >= cs for c, cs, ce in coords) else "other"
        r["set"] = "new"
        r["agarico"] = True
        r["independent"] = False
        rows.append(r)


def line(name, sel):
    k = sum(rule(r) for r in sel)
    n = len(sel)
    lo, hi = wilson(k, n)
    return f"{name}\t{k}/{n}\t{100*k/max(1,n):.1f}%\t[{100*lo:.1f}, {100*hi:.1f}]"


print("group\tflagged/n\trate\tWilson95% (T strict, +-10 kb)")
for name, sel in [
    ("mating, independent (all lineages)", [r for r in rows if r["mating"] == "mating" and r["independent"]]),
    ("mating, independent, Agaricomycetes", [r for r in rows if r["mating"] == "mating" and r["independent"] and r["agarico"]]),
    ("mating, CAAX-selected records (not independent)", [r for r in rows if r["mating"] == "mating" and not r["independent"]]),
    ("other, all lineages", [r for r in rows if r["mating"] == "other"]),
    ("other, Agaricomycetes", [r for r in rows if r["mating"] == "other" and r["agarico"]]),
    ("other, old panel", [r for r in rows if r["mating"] == "other" and r["set"] == "old"]),
    ("other, new genomes", [r for r in rows if r["mating"] == "other" and r["set"] == "new"]),
]:
    print(line(name, sel))

print("\nper genome\tmating\tother\tmating_flagged\tother_flagged")
for asm in sorted({r["asm"] for r in rows}):
    g = [r for r in rows if r["asm"] == asm]
    m = [r for r in g if r["mating"] == "mating"]
    o = [r for r in g if r["mating"] == "other"]
    print(f"{asm}\t{len(m)}\t{len(o)}\t{sum(rule(r) for r in m)}\t{sum(rule(r) for r in o)}")

with open(os.path.join(HERE, "panel_loci.tsv"), "w") as fo:
    keys = ["set", "asm", "contig", "start", "end", "mating", "independent", "agarico", "T_10kb", "T2_10kb", "T_20kb", "Hx_10kb"]
    fo.write("\t".join(keys) + "\n")
    for r in rows:
        fo.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")
