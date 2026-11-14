#!/usr/bin/env python3
"""Infer library strandedness from a BAM and a GFF3, and report read orientation over named regions.

Read-only analysis (needs samtools on PATH). For up to --n annotated genes of at least --min-len bp (and not overlapping another
gene), count properly paired reads by the orientation of read 1 relative to the gene strand:
  same     read 1 on the gene strand  (forward-stranded library, 'fr-secondstrand')
  opposite read 1 on the opposite strand (reverse-stranded library, dUTP, 'fr-firststrand')
A library is called stranded when one orientation holds >= --threshold of the reads; otherwise unstranded.
For each region (name:contig:start:end) the same count is reported as read-1 forward / reverse counts and, if the library is
stranded, the implied transcript strand.
Usage: infer_strandedness.py --bam X.bam --gff3 X.gff3 [--region name:contig:start:end ...]
"""
import argparse
import gzip
import subprocess


def classify(same, opposite, threshold=0.85):
    """-> 'forward-stranded', 'reverse-stranded' or 'unstranded' from read counts relative to gene strand."""
    tot = same + opposite
    if tot == 0:
        return "no_reads"
    if same / tot >= threshold:
        return "forward-stranded"
    if opposite / tot >= threshold:
        return "reverse-stranded"
    return "unstranded"


def implied_strand(r1_fwd, r1_rev, library):
    """Transcript strand implied by read-1 orientation in a stranded library."""
    if library == "forward-stranded":
        return "+" if r1_fwd > r1_rev else "-"
    if library == "reverse-stranded":
        return "-" if r1_fwd > r1_rev else "+"
    return "unknown"


def reads_r1(bam, region):
    """(r1_forward, r1_reverse) among primary, mapped, properly paired reads in a region."""
    out = subprocess.run(["samtools", "view", "-f", "66", "-F", "3852", bam, region], capture_output=True, text=True).stdout
    fw = rv = 0
    for line in out.splitlines():
        flag = int(line.split("\t", 2)[1])
        if flag & 64:
            if flag & 16:
                rv += 1
            else:
                fw += 1
    return fw, rv


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--region", action="append", default=[])
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--min-len", type=int, default=1000)
    ap.add_argument("--threshold", type=float, default=0.85)
    a = ap.parse_args()
    genes = []
    op = gzip.open if a.gff3.endswith(".gz") else open
    with op(a.gff3, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] == "gene":
                genes.append((f[0], int(f[3]), int(f[4]), f[6]))
    genes.sort()
    keep = []
    for i, g in enumerate(genes):
        prev_end = genes[i - 1][2] if i and genes[i - 1][0] == g[0] else -1
        next_start = genes[i + 1][1] if i + 1 < len(genes) and genes[i + 1][0] == g[0] else 10 ** 12
        if g[2] - g[1] + 1 >= a.min_len and g[1] > prev_end and g[2] < next_start:
            keep.append(g)
    step = max(1, len(keep) // a.n)
    same = opp = 0
    for c, s, e, st in keep[::step][:a.n]:
        fw, rv = reads_r1(a.bam, f"{c}:{s}-{e}")
        same += fw if st == "+" else rv
        opp += rv if st == "+" else fw
    lib = classify(same, opp, a.threshold)
    print(f"genes used {len(keep[::step][:a.n])}; read1 same-strand {same}, opposite {opp}; share same {same / max(1, same + opp):.3f}; library: {lib}")
    for r in a.region:
        name, c, s, e = r.split(":")
        fw, rv = reads_r1(a.bam, f"{c}:{s}-{e}")
        print(f"{name}\tread1_forward {fw}\tread1_reverse {rv}\timplied_strand {implied_strand(fw, rv, lib)}")


if __name__ == "__main__":
    main()
