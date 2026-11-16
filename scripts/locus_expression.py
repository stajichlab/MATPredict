#!/usr/bin/env python3
"""Expression of named loci relative to all annotated genes, from `samtools bedcov` output.

Read-only analysis. Inputs
  --gff3      BFD GFF3 (gene features)
  --bedcov    output of `samtools bedcov REGIONS.bed BAM` (chrom start end name sum_of_depth) for every gene and every locus region
  --loci      TSV name<TAB>chrom<TAB>start<TAB>end (1-based inclusive) of the extra regions that were added to the BED, plus
              --genes ids of annotated genes to report
Output (stdout and --out): for each reported locus: length, mean depth, depth percentile among all annotated genes with length >= 300 bp,
and the share of annotated genes with mean depth >= 1.
Mean depth is not TPM: it is not corrected for library size; compare within one library only.
"""
import argparse


def parse_bedcov(path):
    """-> {name: (length, mean_depth)} with BED half-open coordinates."""
    out = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        s, e = int(f[1]), int(f[2])
        out[f[3]] = (e - s, float(f[4]) / max(1, e - s))
    return out


def percentile_of(value, values):
    """Share of values strictly below `value`."""
    return sum(v < value for v in values) / len(values) if values else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bedcov", required=True)
    ap.add_argument("--report", nargs="+", required=True, help="names in the bedcov file to report")
    ap.add_argument("--gene-prefix", default="", help="names that start with this prefix are annotated genes (the background)")
    ap.add_argument("--min-len", type=int, default=300)
    ap.add_argument("--out")
    a = ap.parse_args()
    cov = parse_bedcov(a.bedcov)
    bg = [d for n, (l, d) in cov.items() if n.startswith(a.gene_prefix) and l >= a.min_len and n not in a.report]
    lines = ["name\tlength\tmean_depth\tpercentile_among_genes\tgenes_with_depth_ge_1"]
    for n in a.report:
        if n not in cov:
            lines.append(f"{n}\tNA\tNA\tNA\tNA")
            continue
        l, d = cov[n]
        lines.append(f"{n}\t{l}\t{d:.2f}\t{percentile_of(d, bg):.3f}\t{sum(x >= 1 for x in bg) / len(bg):.3f}")
    text = "\n".join(lines) + "\n"
    print(text)
    if a.out:
        open(a.out, "w").write(text)


if __name__ == "__main__":
    main()
