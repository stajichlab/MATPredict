#!/usr/bin/env python
"""Stratified random sample of v0.6.0 Ascomycota genomes with status `called` (they have at least one locus report).

Strata = (class_, order) of genomes.tsv (114). n_h = min(N_h, max(FLOOR, round(f * N_h))), f chosen so sum n_h ~ TARGET. Seed fixed.
Writes sample.tsv (genome, class_, order, N_h, n_h, weight = N_h / n_h) for extrapolation.
"""
import csv, random, sys
from collections import defaultdict
TARGET, FLOOR, SEED = 1200, 5, 20261017
rows = [r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t") if r["status"] == "called"]
strata = defaultdict(list)
for r in rows:
    strata[(r["class_"], r["order"])].append(r["genome"])
N = len(rows)
def alloc(f):
    return {k: min(len(v), max(FLOOR, round(f * len(v)))) for k, v in strata.items()}
lo, hi = 0.0, 1.0
for _ in range(50):
    f = (lo + hi) / 2
    if sum(alloc(f).values()) > TARGET: hi = f
    else: lo = f
a = alloc(lo)
rng = random.Random(SEED)
with open(sys.argv[2], "w") as fh:
    fh.write("genome\tclass_\torder\tN_h\tn_h\tweight\n")
    for k in sorted(strata):
        for g in sorted(rng.sample(sorted(strata[k]), a[k])):
            fh.write(f"{g}\t{k[0]}\t{k[1]}\t{len(strata[k])}\t{a[k]}\t{len(strata[k]) / a[k]:.4f}\n")
print(N, "called genomes;", sum(a.values()), "sampled;", len(strata), "strata; f =", round(lo, 4))
