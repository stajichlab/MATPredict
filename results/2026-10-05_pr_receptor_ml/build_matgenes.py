#!/usr/bin/env python3
"""Collect the non-receptor, non-precursor MAT gene proteins of the Basidiomycota
records (HD-class genes and conserved flank genes) into matgenes.faa.
Header: >record|class|name   class = HD or FLANK.
Usage: build_matgenes.py REPO_ROOT OUT.faa"""
import glob
import sys

ROOT = sys.argv[1]
HD = {"HD1", "HD2", "bE", "bW", "SXI1", "SXI2"}
FLANK = {"MIP1", "beta_fg", "STE20", "UAP1", "FAO1", "RPO41", "NOG2", "PAN6", "RPL39", "SLA2"}
out = open(sys.argv[2], "w")
n = {"HD": 0, "FLANK": 0}
for f in sorted(glob.glob(ROOT + "/db/Basidiomycota/*/*/proteins.faa")):
    rec = f.split("/")[-2]
    cls = None
    for line in open(f):
        if line.startswith(">"):
            name = None
            for tok in line[1:].split("|"):
                if tok.startswith("name="):
                    name = tok[5:].strip()
            cls = "HD" if name in HD else ("FLANK" if name in FLANK else None)
            if cls:
                out.write(f">{rec}|{cls}|{name}\n")
                n[cls] += 1
        elif cls:
            out.write(line)
print(n)
