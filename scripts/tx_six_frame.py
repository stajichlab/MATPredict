#!/usr/bin/env python3
"""Six-frame translation of the sequences in a nucleotide FASTA (transcripts), one protein record per frame.
Read-only helper. Header: '<id>|sense|f1' ... '<id>|anti|f3'. Stop codons are '*'; non-ACGT codons give 'X'.
Usage: tx_six_frame.py --fasta TX.fa[.gz] --out OUT.faa
"""
import argparse
import gzip
import sys

BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(seq):
    return "".join(CODON.get(seq[i:i + 3].upper(), "X") for i in range(0, len(seq) - 2, 3))


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    with open(a.out, "w") as fo:
        for name, seq in read_fasta(a.fasta):
            for label, s in (("sense", seq), ("anti", revcomp(seq))):
                for f in range(3):
                    fo.write(f">{name}|{label}|f{f + 1}\n{translate(s[f:])}\n")


if __name__ == "__main__":
    sys.exit(main())
