#!/usr/bin/env python3
"""Per-window depth of a panel FASTA's k-mers in one strain's reads (first N reads): where does an idiomorph panel pick up depth?
usage: profile_panel.py panel.fasta reads_1.fq.gz reads_2.fq.gz [max_reads] [window]"""
import gzip, sys
K = 31
def rc(s): return s[::-1].translate(str.maketrans("ACGT", "TGCA"))
seq = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">")).upper()
maxr = int(sys.argv[4]) if len(sys.argv) > 4 else 4_000_000; win = int(sys.argv[5]) if len(sys.argv) > 5 else 250
pos = {}
for i in range(len(seq) - K + 1):
    k = seq[i:i+K]; pos.setdefault(min(k, rc(k)), []).append(i)
cnt = [0] * (len(seq) - K + 1)
for f in sys.argv[2:4]:
    with gzip.open(f, "rt") as fh:
        for n, line in enumerate(fh):
            if n % 4 == 1:
                if n // 4 >= maxr: break
                r = line.strip().upper()
                for j in range(len(r) - K + 1):
                    k = r[j:j+K]
                    for p in pos.get(min(k, rc(k)), ()): cnt[p] += 1
for s in range(0, len(cnt), win):
    seg = cnt[s:s+win]; print(s, s + len(seg), round(sum(seg) / len(seg), 1))
