#!/usr/bin/env python3
"""Write the proteins of an Augustus GFF3 (run with --protein=on) as FASTA, one per transcript.

Augustus prints each protein in comment lines of its GFF3:
  # start gene g1
  ...
  # protein sequence = [MKV...
  # ...LLA]
  # end gene g1
The header is rewritten to <gene>-T<n> (g1.t1 -> g1-T1) so flank_first_scan.py maps a protein to its gene.
Usage: augustus_proteins.py IN.gff3 OUT.faa[.gz]
"""
import gzip
import re
import sys


def proteins(lines):
    """Yield (transcript_id, sequence) for every protein block."""
    gene, buf, on = None, [], False
    for line in lines:
        if line.startswith("# start gene "):
            gene = line.split()[-1]
        if "# protein sequence = [" in line:
            on = True
            buf = [line.split("[", 1)[1].strip()]
        elif on:
            buf.append(line[1:].strip())
        if on and "]" in buf[-1]:
            seq = "".join(buf).replace("]", "").replace(" ", "")
            yield gene, seq
            on, buf = False, []


def main():
    op = gzip.open if sys.argv[2].endswith(".gz") else open
    n = 0
    with op(sys.argv[2], "wt") as fo:
        for gene, seq in proteins(open(sys.argv[1])):
            fo.write(f">{gene}-T1\n{seq}\n")
            n += 1
    print(n, "proteins")


if __name__ == "__main__":
    main()
