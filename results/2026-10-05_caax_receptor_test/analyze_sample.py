#!/usr/bin/env python3
"""Observed vs chance strict-CAAX flags on STE3-like loci, per order.

Per genome g (from out_sample/<asm>.loci.tsv and .random.tsv):
  n_g      STE3-like loci found by miniprot
  o_g      loci with a strict-CAAX ORF (T) within 10 kb
  p_g      fraction of 1,000 random windows (same size, away from STE3 loci) with a T ORF
  e_g      n_g * p_g, the number of loci expected to be flagged by chance
Per order: O = sum o_g, E = sum e_g, ratio O/E, and excess fraction (O-E)/O
(the share of flagged loci that exceeds chance). Intervals: bootstrap over genomes.
This is a population-level estimate; it does not say which single locus is real.
"""
import csv
import os
import random
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(sys.argv[0]))
SAMPLE = os.path.join(HERE, "sample_genomes.tsv")
OUT = os.path.join(HERE, "out_sample")
N_BOOT = 2000
SEED = 7


def read(f):
    return list(csv.DictReader(open(f), delimiter="\t"))


per_genome = []
missing = []
for s in read(SAMPLE):
    a = s["genome"]
    lf, rf = os.path.join(OUT, a + ".loci.tsv"), os.path.join(OUT, a + ".random.tsv")
    if not (os.path.exists(lf) and os.path.exists(rf)):
        missing.append(a)
        continue
    loci, rnd = read(lf), read(rf)
    n = len(loci)
    o = sum(int(r["T_10kb"]) > 0 for r in loci)
    p = sum(int(r["T_10kb"]) > 0 for r in rnd) / max(1, len(rnd))
    per_genome.append({
        "genome": a, "order": s["order"] or "(no order)", "subphylum": s["subphylum"],
        "has_unverified_call": int(s["has_unverified_call"]),
        "n_loci": n, "flagged": o, "p_random": p, "expected": n * p,
    })

with open(os.path.join(HERE, "sample_per_genome.tsv"), "w") as fo:
    keys = list(per_genome[0])
    fo.write("\t".join(keys) + "\n")
    for r in per_genome:
        fo.write("\t".join(f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")


def stat(gs, rng=None):
    if rng:
        gs = [rng.choice(gs) for _ in gs]
    O = sum(g["flagged"] for g in gs)
    E = sum(g["expected"] for g in gs)
    return O, E, (O / E if E else float("nan")), ((O - E) / O if O else float("nan"))


def boot(gs):
    rng = random.Random(SEED)
    xs = sorted(x for x in (stat(gs, rng)[3] for _ in range(N_BOOT)) if x == x)
    if not xs:
        return (float("nan"), float("nan"))
    return xs[int(0.025 * len(xs))], xs[int(0.975 * len(xs)) - 1]


by_order = defaultdict(list)
for g in per_genome:
    by_order[g["order"]].append(g)

lines = ["order\tgenomes\tloci\tflagged_O\texpected_E\tO/E\texcess_(O-E)/O\tboot95_lo\tboot95_hi"]
for order, gs in sorted(by_order.items(), key=lambda kv: -sum(g["flagged"] for g in kv[1])):
    O, E, ratio, ex = stat(gs)
    lo, hi = boot(gs)
    lines.append(f"{order}\t{len(gs)}\t{sum(g['n_loci'] for g in gs)}\t{O}\t{E:.1f}\t{ratio:.2f}\t{ex:.2f}\t{lo:.2f}\t{hi:.2f}")
O, E, ratio, ex = stat(per_genome)
lo, hi = boot(per_genome)
lines.append(f"ALL\t{len(per_genome)}\t{sum(g['n_loci'] for g in per_genome)}\t{O}\t{E:.1f}\t{ratio:.2f}\t{ex:.2f}\t{lo:.2f}\t{hi:.2f}")

# genome-level: does a flagged locus predict the v0.6.0 unverified call?
tab = defaultdict(int)
for g in per_genome:
    tab[(g["flagged"] > 0, bool(g["has_unverified_call"]))] += 1
lines.append("")
lines.append("genomes: any flagged locus (rows) vs unverified PR call in v0.6.0 (columns)")
lines.append("flagged\tcalled\tuncalled")
for fl in (True, False):
    lines.append(f"{fl}\t{tab[(fl, True)]}\t{tab[(fl, False)]}")
# Apply each order's chance share (E/O, capped at 1) to its unverified calls.
# Approximation: calls are merged loci, flagged loci are miniprot loci.
CALLS = os.path.join(HERE, "unverified_calls_by_order.tsv")
if os.path.exists(CALLS):
    lines.append("")
    lines.append("1,360 unverified calls by order, with the chance-expected number")
    lines.append("order\tcalls\tchance_share\tchance_calls\tbeyond_chance_calls\tgenomes_sampled")
    tot_calls = tot_chance = 0.0
    for row in read(CALLS):
        order = row["order"] or "(no order)"
        n_calls = int(row["calls"])
        gs = by_order.get(order, [])
        O, E, _, _ = stat(gs) if gs else (0, 0, 0, 0)
        share = min(1.0, E / O) if O else float("nan")
        chance = n_calls * share if share == share else float("nan")
        lines.append(f"{order}\t{n_calls}\t{share:.2f}\t{chance:.0f}\t{n_calls - chance:.0f}\t{len(gs)}")
        if chance == chance:
            tot_calls += n_calls
            tot_chance += chance
    lines.append(f"TOTAL\t{tot_calls:.0f}\t{tot_chance / tot_calls:.2f}\t{tot_chance:.0f}\t{tot_calls - tot_chance:.0f}\t")
if missing:
    lines.append(f"\nMISSING scans: {len(missing)}: " + ",".join(missing[:10]))

text = "\n".join(lines)
print(text)
open(os.path.join(HERE, "sample_summary.txt"), "w").write(text + "\n")
