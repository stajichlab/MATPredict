#!/usr/bin/env python3
"""reads-type panel for a HOMOTHALLIC species whose MAT genes sit in ONE locus (Sclerotinia sclerotiorum): two 'idiomorphs' cut from the same locus.

The MAT1-1 file holds the core of the alpha-box gene (MAT1-1-1) and the MAT1-2 file the core of the HMG gene (MAT1-2-1); the flanks and the region between the genes
go into both files, so they are the shared control and only the two cores are specific. Every strain of a homothallic species should therefore be called `both`
(a positive control of the typing, not a mating-type result); a one-idiomorph call would mean one of the genes is missing or too divergent.

    panel_homothallic.py --genome-fa G.fa.gz --contig C --span START-END --core1 START-END --core2 START-END --out DIR --prefix sclerotinia [--margin 300 --max-flank 6000]
Coordinates are 1-based, inclusive, on the contig (the `detect` locus span and the gene extents from its report).
"""
import argparse, gzip, os


def contig(path, name):
    seq = []; on = False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if on: break
                on = line[1:].split()[0] == name; continue
            if on: seq.append(line.strip())
    return "".join(seq)


def rng(s):
    a, b = s.split("-"); return int(a), int(b)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for k in ("genome-fa", "contig", "span", "core1", "core2", "out", "prefix"): ap.add_argument("--" + k, required=True)
    ap.add_argument("--margin", type=int, default=300); ap.add_argument("--max-flank", type=int, default=6000)
    a = ap.parse_args()
    sq = contig(a.genome_fa, a.contig); s0, s1 = rng(a.span); c1 = rng(a.core1); c2 = rng(a.core2); m = a.margin
    cut = lambda x, y: sq[x - 1:y]
    c1 = (c1[0] - m, c1[1] + m); c2 = (c2[0] - m, c2[1] + m)
    lo, hi = sorted([c1, c2])
    left = cut(max(s0, lo[0] - a.max_flank), lo[0] - 1); mid = cut(lo[1] + 1, hi[0] - 1); right = cut(hi[1] + 1, min(s1, hi[1] + a.max_flank))
    os.makedirs(a.out, exist_ok=True)
    def wr(fh, name, seq):
        if len(seq) >= 200: fh.write(f">{name}\n" + "\n".join(seq[j:j + 80] for j in range(0, len(seq), 80)) + "\n")
    for idio, core in (("MAT1-1", c1), ("MAT1-2", c2)):
        with open(os.path.join(a.out, f"{a.prefix}_{idio}.fasta"), "w") as fh:
            wr(fh, f"{a.prefix}_{idio}_core {a.contig}:{core[0]}-{core[1]}", cut(*core))
            wr(fh, f"{a.prefix}_flank_L {a.contig}", left); wr(fh, f"{a.prefix}_flank_M {a.contig}", mid); wr(fh, f"{a.prefix}_flank_R {a.contig}", right)
    print("core1", c1, "core2", c2, "flank lengths", len(left), len(mid), len(right))


if __name__ == "__main__":
    main()
