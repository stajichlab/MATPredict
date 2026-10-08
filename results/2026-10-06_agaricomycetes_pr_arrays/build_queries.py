#!/usr/bin/env python3
"""Region queries from curated db records (Basidiomycota): class HD (HD1/HD2/bE/bW),
STE20, MIPBF (MIP1, beta_fg: flank genes of the A/HD locus in Agaricomycetes).
Header >record|class|name. Usage: build_queries.py REPO_ROOT OUT.faa"""
import glob, sys
ROOT, OUT = sys.argv[1], sys.argv[2]
CLS = {"HD1": "HD", "HD2": "HD", "bE": "HD", "bW": "HD", "STE20": "STE20", "MIP1": "MIPBF", "beta_fg": "MIPBF"}
AGARI = {"Agaricales", "Polyporales", "Russulales", "Boletales"}   # HD and MIP1/beta_fg queries: Agaricomycete records only
n = {}
with open(OUT, "w") as fo:
    for f in sorted(glob.glob(ROOT + "/db/Basidiomycota/*/*/proteins.faa")):
        rec, cls = f.split("/")[-2], None
        for line in open(f):
            if line.startswith(">"):
                name = next((t[5:].strip() for t in line[1:].split("|") if t.startswith("name=")), None)
                cls = CLS.get(name)
                if cls in ("HD", "MIPBF") and f.split("/")[-3] not in AGARI:
                    cls = None
                if cls:
                    fo.write(f">{rec}|{cls}|{name}\n"); n[cls] = n.get(cls, 0) + 1
            elif cls:
                fo.write(line)
print(n)
