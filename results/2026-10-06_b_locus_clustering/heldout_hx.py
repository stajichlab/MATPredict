#!/usr/bin/env python3
"""Held-out tblastn precursor homology for the genomes that carry a curated B record.
Reuses scan_genome.homology_hits (same settings: -seg no, word 2, E<=1) and
pheromones_curated.faa, with the records of the genome's own species removed
(pheromones_curated.faa holds Schizophyllum 5334_ and Coprinopsis 5346_ Agaricomycete
precursors; all other curated Agaricomycete B records are not in it).
Usage: heldout_hx.py ASMID EXCLUDE_PREFIX_OR_NONE OUTDIR  (run from the array-study dir)"""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../2026-10-06_agaricomycetes_pr_arrays")
import scan_genome as sg
asm, excl, out = sys.argv[1:4]
src = sg.PHERO
tmp = tempfile.mkdtemp(prefix="hx_", dir=os.environ.get("SCRATCH", "/tmp"))
filt = os.path.join(tmp, "phero_filt.faa")
keep = True
with open(src) as fi, open(filt, "w") as fo:
    for l in fi:
        if l.startswith(">"):
            keep = not (excl != "none" and l[1:].startswith(excl))
        if keep:
            fo.write(l)
sg.PHERO = filt
seqs = sg.read_fasta(os.path.join(sg.LIB, asm + ".fa.gz"))
fa = os.path.join(tmp, asm + ".fa")
with open(fa, "w") as fo:
    for k, v in seqs.items():
        fo.write(f">{k}\n{v}\n")
hits = sg.homology_hits(fa, tmp)
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, asm + ".hx_heldout.tsv"), "w") as fo:
    fo.write("contig\tpos\tquery\tevalue\texcluded\n")
    for c, p, q, e in hits:
        fo.write(f"{c}\t{p}\t{q}\t{e:g}\t{excl}\n")
