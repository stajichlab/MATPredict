#!/usr/bin/env python3
"""Contamination check: protein identity of every gene on the contig that carries a candidate locus
against a reference proteome from another genus. Usage: contig_origin_check.py STAGE1_DIR GENOME GENE REFGENOME
Prints contig length, gene count, and the identity distribution of the contig's genes vs the reference."""
import os, statistics, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from homeodomain_tree import Prot, table
from flank_first_scan import parse_gff
ANN = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation"
d, g, gene, ref = sys.argv[1:5]
P = Prot()
rows = table(f"{d}/{g}.genes_all.tsv")
ctg = next(r[1] for r in rows if r[0] == gene)
genes = [r for r in rows if r[1] == ctg]
print(g, "contig", ctg, "genes", len(genes), "span", min(int(r[2]) for r in genes), max(int(r[3]) for r in genes))
t = os.environ.get("SCRATCH", "/tmp")
q = os.path.join(t, "q.faa"); r_ = os.path.join(t, "r.faa")
with open(q, "w") as fo:
    for r in genes:
        s = P.get(g).get(r[0])
        if s: fo.write(f">{r[0]}\n{s}\n")
with open(r_, "w") as fo:
    for k, s in P.get(ref).items(): fo.write(f">{k}\n{s}\n")
subprocess.run(["diamond", "makedb", "--in", r_, "-d", r_, "--quiet"], check=True)
subprocess.run(["diamond", "blastp", "-q", q, "-d", r_ + ".dmnd", "--max-target-seqs", "1", "--evalue", "1e-10",
                "--outfmt", "6", "qseqid", "pident", "qcovhsp", "-o", q + ".tsv", "--quiet"], check=True)
ids = [float(l.split("\t")[1]) for l in open(q + ".tsv")]
print("genes with a hit vs", ref, ":", len(ids), "of", len(genes))
if ids:
    print("identity median %.1f  min %.1f  max %.1f; genes >=95%%: %d" % (statistics.median(ids), min(ids), max(ids), sum(i >= 95 for i in ids)))
