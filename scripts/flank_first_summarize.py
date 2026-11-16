#!/usr/bin/env python3
"""Summarise stage-1 loci per lineage: which genomes carry MAT-like architectures.

Architectures counted per genome (a locus = candidates within WINDOW genes):
  HD1+HD2     locus with control hits of class HD1 and HD2 (Ropars 2016 AMF proteins)
  sexM/sexP   locus with a control hit of class sexM or sexP (Mucorales HMG)
  HMG+flank   sexM/sexP locus that also has a tptA/rnhA/glrA/algA/btbA control hit
  2xHomeo     locus with two or more genes with a Pfam Homeodomain hit
  HMG+Homeo   locus with a Pfam HMG_box (or HMG_box_2) gene and a Pfam Homeodomain gene
Usage: flank_first_summarize.py STAGE1_DIR LIST[:LABEL] [LIST[:LABEL] ...] > table.tsv
"""
import collections, os, sys

d = sys.argv[1]
FL = {"tptA", "rnhA", "glrA", "algA", "btbA"}
print("lineage\tgenomes\twith_stage1\tHD1+HD2\tsexM/sexP\tHMG+flank\t2xHomeo\tHMG+Homeo")
rows = {}
for spec in sys.argv[2:]:
    path, _, label = spec.partition(":")
    label = label or os.path.basename(path)
    gs = [x.strip() for x in open(path) if x.strip()]
    tot = collections.Counter()
    for g in gs:
        lf, cf = f"{d}/{g}.loci.tsv", f"{d}/{g}.candidates.tsv"
        if not os.path.exists(lf):
            continue
        tot["done"] += 1
        cand = [l.rstrip("\n").split("\t") for l in open(cf)][1:]
        loci = [l.rstrip("\n").split("\t") for l in open(lf)][1:]
        feat = collections.defaultdict(set)
        for l in loci:
            a, b, ctg = int(l[2]), int(l[3]), l[1]
            hs = [c for c in cand if c[1] == ctg and a <= int(c[2]) <= b]
            cls = {c[5] for c in hs if c[5]}
            pf = [c[3] for c in hs if c[3]]
            if {"HD1", "HD2"} <= cls: feat["HD1+HD2"].add(l[0])
            if {"sexM", "sexP"} & cls: feat["sexM/sexP"].add(l[0])
            if {"sexM", "sexP"} & cls and FL & cls: feat["HMG+flank"].add(l[0])
            if pf.count("Homeodomain") >= 2: feat["2xHomeo"].add(l[0])
            if any(p.startswith("HMG_box") for p in pf) and "Homeodomain" in pf: feat["HMG+Homeo"].add(l[0])
        for k in feat: tot[k] += 1
    print("\t".join(map(str, [label, len(gs), tot["done"]] + [tot[k] for k in ["HD1+HD2", "sexM/sexP", "HMG+flank", "2xHomeo", "HMG+Homeo"]])))
