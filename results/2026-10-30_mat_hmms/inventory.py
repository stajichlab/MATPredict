"""Curated MAT proteins in db/ (Saccharomycotina only) by gene name -> class a1 / a2 / alpha1 / alpha2."""
import re, glob, collections, sys
ROOT = "../../db/Ascomycota"
def cls(n):
    n = n.lower().replace("mtl", "mat")
    m = re.fullmatch(r"mat(a|alpha)([12])", n)
    return (m.group(1) + m.group(2)) if m else None
out = []; cnt = collections.Counter()
for f in sorted(glob.glob(ROOT + "/*/*/proteins.faa")):
    parts = f.split("/")
    order, rec = parts[-3], parts[-2]
    if order in ("Pezizales", "Eurotiales", "Hypocreales", "Helotiales", "Sordariales", "Onygenales", "Glomerellales", "Diaporthales", "Mycosphaerellales", "Microascales", "Lecanorales", "Ophiostomatales", "Xylariales", "Capnodiales", "Pleosporales", "Chaetothyriales", "Sordariomycetes"): continue
    name = None; seq = []
    for l in open(f):
        if l.startswith(">"):
            if name: out.append((rec, name, "".join(seq)))
            m = re.search(r"name=([^|]+)", l); name = m.group(1) if m else l[1:].split("|")[-1].strip(); seq = []
        else: seq.append(l.strip())
    if name: out.append((rec, name, "".join(seq)))
keep = []
for rec, name, s in out:
    c = cls(name)
    if c: keep.append((rec, name, c, s)); cnt[(rec.split("_")[0], c)] += 1
with open("data/curated_mat_proteins.faa", "w") as fh:
    for rec, name, c, s in keep: fh.write(f">{rec}|{name}|{c}\n{s}\n")
orders = collections.Counter()
for rec, name, c, s in keep: orders[c] += 1
print(orders); print(sorted({r.split('_')[0] for r, *_ in keep}))
for rec, name, c, s in keep:
    if c in ("a1", "alpha2", "a2"): print(c, rec, name, len(s))
