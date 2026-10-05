#!/usr/bin/env python3
"""Does homology to curated pheromone precursors (Hx) add to the strict-CAAX rule?

Rules per STE3-like locus, window +-10 kb:
  T      strict-CAAX ORF
  Hx     tblastn hit (E<=1) of the 31 curated precursors (any ORF)
  T&Hx   both
Panel: mating (independent records) vs other, Wilson 95%.
Sample (502 genomes): observed vs chance (random windows of the same size) per
rule, pooled and per order. Genomes of genera that supplied curated precursors
(Cryptococcus, Schizophyllum, Coprinopsis, Ustilaginales genera) are left out of
the sample statistics: their own precursors would match (Hx not held out).
Usage: analyze_hx.py GENOMES_TSV  (v0.6.0 genomes.tsv, for species names)
"""
import csv
import math
import os
import random
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(sys.argv[0]))
GENOMES = sys.argv[1]
SOURCE_GENERA = {"Cryptococcus", "Schizophyllum", "Coprinopsis", "Ustilago", "Mycosarcoma",
                 "Pseudozyma", "Moesziomyces", "Tranzscheliella", "Farysia", "Sporisorium",
                 "Malassezia", "Cryptococcus", "Kwoniella", "Filobasidiella"}
RULES = {"T": lambda t, h: t, "Hx": lambda t, h: h, "T&Hx": lambda t, h: t and h}


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def read(f):
    return list(csv.DictReader(open(f), delimiter="\t"))


out = []
# ---- panel
panel = read(os.path.join(HERE, "panel_loci.tsv"))
out.append("PANEL (+-10 kb)\ngroup\trule\tflagged/n\trate\tWilson95%")
groups = {
    "mating, independent": [r for r in panel if r["mating"] == "mating" and r["independent"] == "True"],
    "mating, CAAX-selected": [r for r in panel if r["mating"] == "mating" and r["independent"] == "False"],
    "other": [r for r in panel if r["mating"] == "other"],
}
for gname, rows in groups.items():
    for rname, fn in RULES.items():
        k = sum(bool(fn(int(r["T_10kb"]) > 0, int(r["Hx_10kb"]) > 0)) for r in rows)
        lo, hi = wilson(k, len(rows))
        out.append(f"{gname}\t{rname}\t{k}/{len(rows)}\t{100*k/max(1,len(rows)):.1f}%\t[{100*lo:.1f}, {100*hi:.1f}]")

# ---- sample
species = {g["genome"]: g["species"] for g in read(GENOMES)}
samp = read(os.path.join(HERE, "sample_genomes.tsv"))
per_genome, dropped = [], 0
for s in samp:
    a = s["genome"]
    genus = species.get(a, "").split(" ")[0]
    if genus in SOURCE_GENERA:
        dropped += 1
        continue
    loci = read(os.path.join(HERE, "out_sample", a + ".loci.tsv"))
    rnd = read(os.path.join(HERE, "out_sample", a + ".random.tsv"))
    rec = {"genome": a, "order": s["order"] or "(no order)", "n": len(loci), "rnd_n": len(rnd)}
    for rname, fn in RULES.items():
        rec["O_" + rname] = sum(bool(fn(int(r["T_10kb"]) > 0, int(r["Hx_10kb"]) > 0)) for r in loci)
        p = sum(bool(fn(int(r["T_10kb"]) > 0, int(r["Hx_10kb"]) > 0)) for r in rnd) / max(1, len(rnd))
        rec["E_" + rname] = len(loci) * p
    per_genome.append(rec)


def stat(gs, rule, rng=None):
    if rng:
        gs = [rng.choice(gs) for _ in gs]
    O = sum(g["O_" + rule] for g in gs)
    E = sum(g["E_" + rule] for g in gs)
    return O, E


def excess_ci(gs, rule):
    rng = random.Random(7)
    xs = []
    for _ in range(2000):
        O, E = stat(gs, rule, rng)
        if O:
            xs.append((O - E) / O)
    xs.sort()
    return (xs[int(0.025 * len(xs))], xs[int(0.975 * len(xs)) - 1]) if xs else (float("nan"),) * 2


out.append(f"\nSAMPLE: {len(per_genome)} genomes ({dropped} source-genus genomes left out), "
           f"{sum(g['n'] for g in per_genome)} STE3-like loci")
out.append("rule\tobserved\texpected_by_chance\tratio\texcess_(O-E)/O\tboot95")
for rname in RULES:
    O, E = stat(per_genome, rname)
    lo, hi = excess_ci(per_genome, rname)
    out.append(f"{rname}\t{O}\t{E:.1f}\t{O/E if E else float('nan'):.1f}\t{(O-E)/O if O else float('nan'):.2f}\t{lo:.2f}-{hi:.2f}")

by_order = defaultdict(list)
for g in per_genome:
    by_order[g["order"]].append(g)
out.append("\nper order (observed/expected; genomes)")
out.append("order\tgenomes\tloci\tT\tHx\tT&Hx")
for order, gs in sorted(by_order.items(), key=lambda kv: -sum(g["O_T"] for g in kv[1]))[:20]:
    cells = [f"{stat(gs, r)[0]}/{stat(gs, r)[1]:.1f}" for r in RULES]
    out.append(f"{order}\t{len(gs)}\t{sum(g['n'] for g in gs)}\t" + "\t".join(cells))

text = "\n".join(out)
print(text)
open(os.path.join(HERE, "hx_summary.txt"), "w").write(text + "\n")
