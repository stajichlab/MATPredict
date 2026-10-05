#!/usr/bin/env python3
"""Stratified random sample of Basidiomycota genomes for the chance-model scan.

Strata: every order that holds >=1 unverified (CAAX-only) PR call in the v0.6.0
run. Up to N_PER_ORDER genomes per order, drawn at random from ALL genomes of
that order regardless of call status (so the observed/chance comparison is not
conditioned on a genome having been called). Genomes >= 500 Mb are left out
(scan time), and genomes missing from the BFD library are dropped.
"""
import csv
import os
import random
import sys
from collections import Counter, defaultdict

SRC = sys.argv[1]                      # v0.6.0 results folder (genomes.tsv, loci.tsv)
OUT = sys.argv[2]
LIB = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"
N_PER_ORDER = 30
MAX_BP = 500_000_000
SEED = 20261005

genomes = list(csv.DictReader(open(os.path.join(SRC, "genomes.tsv")), delimiter="\t"))
loci = list(csv.DictReader(open(os.path.join(SRC, "loci.tsv")), delimiter="\t"))
unver = Counter()
called_unver = set()
for r in loci:
    if r["family_called"].endswith(":PR") and "unverified" in r["verification"]:
        unver[r["order"]] += 1
        called_unver.add(r["genome"])

by_order = defaultdict(list)
for g in genomes:
    if g["order"] not in unver or int(g["size_bp"] or 0) >= MAX_BP:
        continue
    if not os.path.exists(os.path.join(LIB, g["genome"] + ".fa.gz")):
        continue
    by_order[g["order"]].append(g)

rng = random.Random(SEED)
rows = []
for order in sorted(by_order):
    pool = by_order[order]
    pick = pool if len(pool) <= N_PER_ORDER else rng.sample(pool, N_PER_ORDER)
    for g in pick:
        rows.append((g["genome"], order, g["subphylum"], len(pool), unver[order],
                     int(g["genome"] in called_unver)))

with open(OUT, "w") as fo:
    fo.write("genome\torder\tsubphylum\torder_pool\torder_unverified_calls\thas_unverified_call\n")
    for r in rows:
        fo.write("\t".join(map(str, r)) + "\n")
print(len(rows), "genomes in", len(by_order), "orders;",
      sum(r[5] for r in rows), "carry an unverified PR call")
