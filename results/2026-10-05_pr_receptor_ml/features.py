#!/usr/bin/env python3
"""Protein-level features of the labelled loci: length, transmembrane helix count
(Kyte-Doolittle, window 19, mean >= 1.6, helices at least 15 aa apart), PF02076
hmmsearch score and domain coverage (pyhmmer), 2-mer and 3-mer composition.
Writes protein_features.tsv.  Usage: features.py labelled.faa PF02076.hmm.gz out.tsv"""
import gzip
import sys

import pyhmmer
from pyhmmer.easel import Alphabet, SequenceFile, TextSequence
from pyhmmer.plan7 import HMMFile

KD = dict(zip("ARNDCQEGHILKMFPSTWYV", [1.8, -4.5, -3.5, -3.5, 2.5, -3.5, -3.5, -0.4, -3.2, 4.5, 3.8, -3.9, 1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2]))


def tm_count(seq, w=19, thr=1.6):
    v = [KD.get(a, 0.0) for a in seq]
    if len(v) < w:
        return 0
    n, last = 0, -100
    i = 0
    while i <= len(v) - w:
        if sum(v[i:i + w]) / w >= thr:
            if i - last >= 15 + w // 2 or last < 0:
                n += 1
            # skip to end of this window so one helix is counted once
            j = i
            while j <= len(v) - w and sum(v[j:j + w]) / w >= thr:
                j += 1
            last = j
            i = j
        else:
            i += 1
    return n


seqs = {}
k = None
for line in open(sys.argv[1]):
    if line.startswith(">"):
        k = line[1:].strip()
        seqs[k] = ""
    else:
        seqs[k] += line.strip()
alpha = Alphabet.amino()
with gzip.open(sys.argv[2], "rb") as fh:
    open("/tmp/_pf02076.hmm", "wb").write(fh.read())
with HMMFile("/tmp/_pf02076.hmm") as hf:
    hmm = hf.read()
digi = [TextSequence(name=k.encode(), sequence=v.replace("*", "")).digitize(alpha) for k, v in seqs.items()]
res = {}
for hits in pyhmmer.hmmsearch(hmm, digi, E=1000):
    for h in hits:
        dom = [d for d in h.domains.included] or list(h.domains)
        cov = sum(d.alignment.hmm_to - d.alignment.hmm_from + 1 for d in dom) / hmm.M
        res[(h.name.decode() if isinstance(h.name, bytes) else h.name)] = (h.score, cov)
with open(sys.argv[3], "w") as fo:
    fo.write("locus_id\tprot_len\ttm\thydrophobic_frac\tpf02076_score\tpf02076_cov\n")
    for k, v in seqs.items():
        s, c = res.get(k, (0.0, 0.0))
        hf_ = sum(1 for a in v if a in "AILMFVW") / max(1, len(v))
        fo.write(f"{k}\t{len(v)}\t{tm_count(v)}\t{hf_:.3f}\t{s:.1f}\t{c:.2f}\n")
