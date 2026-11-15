#!/usr/bin/env python3
"""Count reads carrying the JEL423 or the CLFT044 base at the divergent single-base sites, from `samtools mpileup` text on stdin.

For every pileup line at a site listed in --sites (TSV: contig, 1-based position, JEL423 base) the reads with the reference base (CLFT044:
'.' or ',') and the reads with the JEL423 base are counted; other bases and indels are ignored. Reads are counted once per site.
Output (--out, gzip TSV): contig, pos, n_clft, n_jel.
Summary (stdout) with --min-depth: sites with enough reads, and the share of those sites by JEL fraction f = n_jel / (n_clft + n_jel):
  clft  f <= 0.2   het  0.2 < f < 0.8   jel  f >= 0.8
Usage: samtools mpileup -l SITES.bed -f REF.fa -q 20 -Q 20 -B BAM | bd_site_allele_counts.py --sites SITES.tsv --out OUT.tsv.gz
"""
import argparse
import gzip
import re
import sys

_INDEL = re.compile(r"[+-](\d+)")


def count_bases(bases, jel):
    """(n_ref, n_jel) in an mpileup base string for one site. `jel` is the JEL423 base (upper case)."""
    nref = njel = 0
    i, n = 0, len(bases)
    jel = jel.upper()
    while i < n:
        ch = bases[i]
        if ch == "^":
            i += 2
            continue
        if ch in "$*<>#":
            i += 1
            continue
        if ch in "+-":
            m = _INDEL.match(bases, i)
            i = m.end() + int(m.group(1))
            continue
        if ch in ".,":
            nref += 1
        elif ch.upper() == jel:
            njel += 1
        i += 1
    return nref, njel


def classify(nref, njel, min_depth=5):
    """'clft', 'het', 'jel' or None (too few reads)."""
    tot = nref + njel
    if tot < min_depth:
        return None
    f = njel / tot
    return "clft" if f <= 0.2 else "jel" if f >= 0.8 else "het"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-depth", type=int, default=5)
    a = ap.parse_args()
    jel = {}
    for line in open(a.sites):
        c, p, b = line.rstrip("\n").split("\t")
        jel[(c, int(p))] = b
    tally = {"clft": 0, "het": 0, "jel": 0}
    n_sites = 0
    with gzip.open(a.out, "wt") as fo:
        fo.write("contig\tpos\tn_clft\tn_jel\n")
        for line in sys.stdin:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            key = (f[0], int(f[1]))
            if key not in jel:
                continue
            nref, njel = count_bases(f[4], jel[key])
            fo.write(f"{f[0]}\t{f[1]}\t{nref}\t{njel}\n")
            c = classify(nref, njel, a.min_depth)
            if c:
                tally[c] += 1
                n_sites += 1
    tot = max(1, n_sites)
    print(f"sites with >= {a.min_depth} reads: {n_sites}; clft {tally['clft'] / tot:.3f}  het {tally['het'] / tot:.3f}  jel {tally['jel'] / tot:.3f}")


if __name__ == "__main__":
    main()
