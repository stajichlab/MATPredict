#!/usr/bin/env python3
"""Observed vs chance-expected number of genomes per lineage for adjacency architectures.

Chance model, per genome: n genes of one kind and m genes of another are placed at random among
N genes; two genes are "adjacent" when within WINDOW genes (2*WINDOW neighbours). Expected
P(at least one adjacent pair) = 1 - exp(-n*m*2*W/N) (for the same kind, C(n,2) pairs).
The model ignores contig breaks, which makes it over-estimate chance adjacency. A lineage
whose observed count is near the expected count shows no enrichment.
Usage: flank_first_null.py STAGE1_DIR WINDOW LIST[:LABEL] ...
"""
import math, os, sys, collections
d, W = sys.argv[1], int(sys.argv[2])
FL = {"tptA", "rnhA", "glrA", "algA", "btbA"}
print("lineage\tarchitecture\tgenomes\tobserved\texpected_by_chance")
for spec in sys.argv[3:]:
    path, _, label = spec.partition(":")
    gs = [x.strip() for x in open(path) if x.strip()]
    acc = collections.defaultdict(lambda: [0, 0.0])  # arch -> [obs, exp]
    nd = 0
    for g in gs:
        cf = f"{d}/{g}.candidates.tsv"
        if not os.path.exists(cf): continue
        nd += 1
        N = sum(1 for _ in open(f"{d}/{g}.genes_all.tsv")) - 1
        cand = [l.rstrip("\n").split("\t") for l in open(cf)][1:]
        hmg = {c[0] for c in cand if c[3].startswith("HMG_box")}
        hom = {c[0] for c in cand if c[3] == "Homeodomain"}
        mp = {c[0] for c in cand if c[5] in ("sexM", "sexP")}
        fl = {c[0] for c in cand if c[5] in FL}
        hd1 = {c[0] for c in cand if c[5] == "HD1"}; hd2 = {c[0] for c in cand if c[5] == "HD2"}
        idx = {c[0]: (c[1], int(c[2])) for c in cand}
        def adj(A, B):
            return any(a != b and idx[a][0] == idx[b][0] and abs(idx[a][1] - idx[b][1]) <= W for a in A for b in B)
        tests = {
            "sexM/sexP+flank": (mp, fl, len(mp) * len(fl)),
            "2xHomeodomain": (hom, hom, len(hom) * (len(hom) - 1) / 2),
            "HMG+Homeodomain": (hmg, hom, len(hmg) * len(hom)),
            "HD1+HD2 (Ropars)": (hd1, hd2, len(hd1) * len(hd2)),
        }
        for k, (A, B, pairs) in tests.items():
            acc[k][0] += int(adj(A, B))
            acc[k][1] += 1 - math.exp(-pairs * 2 * W / max(N, 1))
    for k, (o, e) in acc.items():
        print(f"{label or path}\t{k}\t{nd}\t{o}\t{e:.1f}")
