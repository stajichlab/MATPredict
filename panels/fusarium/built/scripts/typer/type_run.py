#!/usr/bin/env python3
"""Standalone driver around the repo's MATPredict.reads (copied unchanged: kmers.py, panel.py, typing.py) that writes the same
columns as `matpredict reads-type`.  usage: type_run.py SAMPLE OUT.tsv NAME=FASTA NAME=FASTA -- R1 R2 [-- --k 31]"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from MATPredict.reads.panel import Panel
from MATPredict.reads.typing import type_reads
sample, out = sys.argv[1], sys.argv[2]
i = sys.argv.index("--")
fastas = dict(a.split("=", 1) for a in sys.argv[3:i])
reads = sys.argv[i+1:]
panel = Panel.from_fastas(fastas, k=31, min_run=100)
res = type_reads(panel, reads)
names = list(fastas)
new = not os.path.exists(out)
with open(out, "a") as fh:
    if new:
        fh.write("\t".join(["sample", "call"] + [f"{n}_breadth" for n in names] + [f"{n}_depth" for n in names] + ["reads_used", "shared_depth", "flags"]) + "\n")
    fh.write("\t".join([sample, res.call] + [f"{res.breadth[n]:.4f}" for n in names] + [f"{res.depth[n]:.2f}" for n in names] + [str(res.reads_used), f"{res.shared_depth:.2f}", ",".join(res.flags) or "-"]) + "\n")
print(sample, res.call, res.breadth, res.depth, res.shared_depth)
