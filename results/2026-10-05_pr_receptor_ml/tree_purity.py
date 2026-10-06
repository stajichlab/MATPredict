#!/usr/bin/env python3
"""Clade placement: does the nearest relative of a locus in an STE3 tree carry the same label?

Reads tree.nwk (FastTree, complete-model loci plus the S. cerevisiae Ste3 outgroup),
labelled_loci.tsv. For each locus, the nearest other-genome locus by patristic distance:
 - same order: label agreement (mating->mating, other->other)
 - different order: same
Compares with the expectation under label shuffling within the candidate pool (1000 shuffles).
Writes tree_purity.tsv.  Usage: tree_purity.py tree.nwk tree_map.tsv
"""
import os
import sys

import numpy as np
import pandas as pd
from Bio import Phylo

HERE = os.path.dirname(os.path.abspath(sys.argv[0]))
L = pd.read_csv(f"{HERE}/labelled_loci.tsv", sep="\t").set_index("locus_id")
L = L[L.label.isin(["mating", "other"])]
tree = Phylo.read(sys.argv[1], "newick")
MAP = dict(l.rstrip("\n").split("\t") for l in open(sys.argv[2]))
for t in tree.get_terminals():
    t.name = MAP.get(t.name, t.name)
tips = [t.name for t in tree.get_terminals() if t.name in L.index]
D = {}
for i, a in enumerate(tips):
    for b in tips[i + 1:]:
        D[(a, b)] = D[(b, a)] = tree.distance(a, b)
rng = np.random.default_rng(1)
rows = []
for scope in ("same_order", "other_order"):
    for lab in ("mating", "other"):
        obs, n = 0, 0
        pools = []
        for a in tips:
            if L.label[a] != lab:
                continue
            cand = [b for b in tips if L.asm[b][:15] != L.asm[a][:15] and ((L.order[b] == L.order[a]) == (scope == "same_order"))]
            if not cand:
                continue
            nn = min(cand, key=lambda b: D[(a, b)])
            obs += L.label[nn] == lab
            n += 1
            pools.append((cand, nn))
        # null: shuffle labels among all tips, recompute agreement
        nulls = []
        labs = np.array([L.label[t] for t in tips])
        for _ in range(1000):
            sh = dict(zip(tips, rng.permutation(labs)))
            c = 0
            m = 0
            for a in tips:
                if L.label[a] != lab:
                    continue
                cand = [b for b in tips if L.asm[b][:15] != L.asm[a][:15] and ((L.order[b] == L.order[a]) == (scope == "same_order"))]
                if not cand:
                    continue
                nn = min(cand, key=lambda b: D[(a, b)])
                c += sh[nn] == lab
                m += 1
            nulls.append(c / max(1, m))
            if _ > 200:
                break
        rows.append(dict(scope=scope, label=lab, n_loci=n, nn_same_label=obs, frac=round(obs / max(1, n), 3),
                         null_mean=round(float(np.mean(nulls)), 3), null_p95=round(float(np.percentile(nulls, 95)), 3)))
pd.DataFrame(rows).to_csv(f"{HERE}/tree_purity.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string())
