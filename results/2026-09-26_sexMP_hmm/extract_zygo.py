"""Proteins annotated (funannotate) inside each Zygo truth locus (+-2 kb)."""
import os
T = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-23_zygo_bar/zygo_truth.tsv"
out = open("zygo_locus_proteins.faa", "w")
for line in open(T):
    org, scaf, s, e, idio, fsa = line.rstrip("\n").split("\t")
    s, e = int(s) - 2000, int(e) + 2000
    d = os.path.dirname(fsa)
    gff = [f for f in os.listdir(d) if f.endswith(".gff3")][0]
    prot = [f for f in os.listdir(d) if f.endswith(".proteins.fa")][0]
    want = set()
    for g in open(os.path.join(d, gff)):
        f = g.split("\t")
        if len(f) > 8 and f[0] == scaf and f[2] == "mRNA" and int(f[4]) >= s and int(f[3]) <= e:
            want.add(f[8].split("ID=")[1].split(";")[0])
    seqs, cur = {}, None
    for l in open(os.path.join(d, prot)):
        if l.startswith(">"):
            cur = l[1:].split()[0]; seqs[cur] = []
        else:
            seqs[cur].append(l.strip())
    n = 0
    for w in sorted(want):
        if w in seqs:
            out.write(f">ZYGO|{org}|{idio}|{w}\n{''.join(seqs[w])}\n"); n += 1
    print(org, idio, len(want), n)
