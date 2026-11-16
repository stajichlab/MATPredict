#!/usr/bin/env python3
"""Six-frame translation of genome regions and a Pfam check (HMG_box, HMG_box_2, Homeodomain).

Read-only discovery tool. Used for regions that a tblastn hit points to but that carry no gene model, to ask whether the
region encodes a homeodomain or HMG-box domain.

Input  --regions TSV: name, genome fasta, contig, start, end (1-based, inclusive)
       --hmm     HMM file with the Pfam models (GA cutoffs)
Method: the region is widened by --pad bp each side; all six frames are translated; stop-free stretches of at least
        --min-aa residues are kept as ORFs (a spliced gene can be split over several ORFs; no intron model is used);
        hmmsearch --cut_ga on the ORFs.
Output --out: region_pfam.tsv (name, region, orf_frame, orf_start, orf_end, aa_len, domain, evalue); orf_start/orf_end are
        1-based and inclusive on the forward strand of the contig
"""
import argparse
import os
import subprocess
import sys

BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(seq):
    return "".join(CODON.get(seq[i:i + 3].upper(), "X") for i in range(0, len(seq) - 2, 3))


def orfs(seq, min_aa):
    """Stop-free stretches of >= min_aa residues in the six frames. Returns (frame, start0, end0, aa).
    Coordinates are 0-based on the forward strand of `seq`; frames 1-3 forward, -1..-3 reverse."""
    out = []
    n = len(seq)
    for strand, s in ((1, seq), (-1, revcomp(seq))):
        for f in range(3):
            prot = translate(s[f:])
            pos = 0
            for chunk in prot.split("*"):
                if len(chunk) >= min_aa:
                    a0 = f + pos * 3
                    a1 = a0 + len(chunk) * 3
                    if strand == 1:
                        out.append((f + 1, a0, a1, chunk))
                    else:
                        out.append((-(f + 1), n - a1, n - a0, chunk))
                pos += len(chunk) + 1
    return out


def read_fasta_region(fa, contig, start, end):
    """Sequence of contig[start-1:end] from a (plain or gz) FASTA, streamed."""
    import gzip
    op = gzip.open if fa.endswith(".gz") else open
    buf, on = [], False
    with op(fa, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if on:
                    break
                on = line[1:].split()[0] == contig
            elif on:
                buf.append(line.strip())
    return "".join(buf)[start - 1:end]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--hmm", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pad", type=int, default=600)
    ap.add_argument("--min-aa", type=int, default=40)
    ap.add_argument("--tmp", default=os.environ.get("SCRATCH", "/tmp"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fo = open(os.path.join(a.out, "region_pfam.tsv"), "w")
    fo.write("name\tregion\torf_frame\torf_start\torf_end\taa_len\tdomain\tevalue\n")
    for line in open(a.regions):
        name, fa, contig, s, e = line.rstrip("\n").split("\t")
        s, e = max(1, int(s) - a.pad), int(e) + a.pad
        seq = read_fasta_region(fa, contig, s, e)
        faa = os.path.join(a.tmp, "orfs.faa")
        meta = {}
        with open(faa, "w") as fh:
            for i, (fr, o0, o1, prot) in enumerate(orfs(seq, a.min_aa)):
                fh.write(f">o{i}\n{prot}\n")
                meta[f"o{i}"] = (fr, s + o0, s + o1 - 1, len(prot))  # 1-based, inclusive
        tbl = os.path.join(a.tmp, "orfs.tbl")
        subprocess.run(["hmmsearch", "--cut_ga", "--domtblout", tbl, "-o", "/dev/null", a.hmm, faa], check=True)
        hits = {}
        for l in open(tbl):
            if l[0] != "#":
                f = l.split()
                hits.setdefault(f[0], []).append((f[3], f[6]))
        region = f"{contig}:{s}-{e}"
        if not hits:
            fo.write(f"{name}\t{region}\t\t\t\t\tnone\t\n")
        for oid, hs in hits.items():
            fr, o0, o1, n = meta[oid]
            for dom, ev in hs:
                fo.write(f"{name}\t{region}\t{fr}\t{o0}\t{o1}\t{n}\t{dom}\t{ev}\n")
    fo.close()
    print(open(os.path.join(a.out, "region_pfam.tsv")).read())


if __name__ == "__main__":
    sys.exit(main())
