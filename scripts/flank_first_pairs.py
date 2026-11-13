#!/usr/bin/env python3
"""Adjacent homeodomain-gene pairs: orientation and spacing, per lineage.

A pair = two genes within WINDOW genes on one contig that are both (a) Pfam Homeodomain hits or
(b) control classes HD1 and HD2. For the closest pair per locus, report strand relation
(divergent -/+ , convergent +/-, same) and the number of genes between them.
Usage: flank_first_pairs.py STAGE1_DIR WINDOW LIST[:LABEL] ...
"""
import collections, os, sys
d, W = sys.argv[1], int(sys.argv[2])
print("lineage\tpair_type\tgenomes_with_pair\tpairs\tdivergent\tconvergent\tsame_strand\tadjacent_genes(0 between)")
for spec in sys.argv[3:]:
    path, _, label = spec.partition(":")
    gs = [x.strip() for x in open(path) if x.strip()]
    acc = {"homeo+homeo": collections.Counter(), "HD1+HD2": collections.Counter()}
    ng = collections.Counter()
    for g in gs:
        cf = f"{d}/{g}.candidates.tsv"
        if not os.path.exists(cf): continue
        strand, start = {}, {}
        for l in list(open(f"{d}/{g}.genes_all.tsv"))[1:]:
            f = l.rstrip("\n").split("\t"); strand[f[0]] = f[4]; start[f[0]] = int(f[2])
        cand = [l.rstrip("\n").split("\t") for l in open(cf)][1:]
        for typ, A, B in [("homeo+homeo", [c for c in cand if c[3] == "Homeodomain"], None),
                          ("HD1+HD2", [c for c in cand if c[5] == "HD1"], [c for c in cand if c[5] == "HD2"])]:
            B = A if B is None else B
            seen = set(); found = False
            for a in A:
                for b in B:
                    if a[0] == b[0] or a[1] != b[1] or abs(int(a[2]) - int(b[2])) > W: continue
                    k = tuple(sorted((a[0], b[0])))
                    if k in seen: continue
                    seen.add(k); found = True
                    l_, r_ = sorted((a, b), key=lambda c: start[c[0]])
                    s1, s2 = strand[l_[0]], strand[r_[0]]
                    rel = "divergent" if (s1, s2) == ("-", "+") else "convergent" if (s1, s2) == ("+", "-") else "same"
                    acc[typ][rel] += 1; acc[typ]["n"] += 1
                    if abs(int(a[2]) - int(b[2])) == 1: acc[typ]["adj"] += 1
            if found: ng[typ] += 1
    for typ, c in acc.items():
        if c["n"]: print(f"{label or path}\t{typ}\t{ng[typ]}\t{c['n']}\t{c['divergent']}\t{c['convergent']}\t{c['same']}\t{c['adj']}")
