#!/usr/bin/env python3
"""6-frame translate FASTA records into ORF segments (>= MINAA aa between stops) for a diamond blastx recruit db.
usage: orfdb.py MINAA out.faa in1.fasta [in2.fasta ...]; record id = <file stem>|<record id>|f<frame>|<nt start>"""
import sys, itertools
from pathlib import Path
b = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
T = {a + b_ + c: aas[i] for i, (a, b_, c) in enumerate(itertools.product(b, b, b))}
def rc(s): return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
def fa(p):
    n, s = None, []
    for l in open(p):
        if l.startswith(">"):
            if n: yield n, "".join(s)
            n, s = l[1:].split()[0], []
        else: s.append(l.strip().upper())
    if n: yield n, "".join(s)
minaa = int(sys.argv[1]); out = open(sys.argv[2], "w"); n = 0
for p in sys.argv[3:]:
    for rid, seq in fa(p):
        for strand, sq in (("+", seq), ("-", rc(seq))):
            for f in range(3):
                prot = "".join(T.get(sq[i:i+3], "X") for i in range(f, len(sq) - 2, 3))
                pos = 0
                for seg in prot.split("*"):
                    if len(seg) >= minaa:
                        out.write(f">{Path(p).stem}|{rid}|{strand}{f}|{f + 3*pos}\n{seg}\n"); n += 1
                    pos += len(seg) + 1
print(n, "ORF segments", file=sys.stderr)
