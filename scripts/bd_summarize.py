#!/usr/bin/env python3
"""Summarise Batrachochytrium read-depth and heterozygosity per isolate.

Inputs (from scripts/bd_coverage_isolates.slurm), in BASE/per_run_<ref>/:
  <run>.regions.bed.gz   mosdepth mean depth in 10 kb windows
  <run>.snps.tsv.gz      CHROM POS REF ALT QUAL DP AD (bcftools call -mv, SNPs, QUAL >= 30)
and BASE/bd_runs.tsv (group, run, isolate, country, ...).

Per isolate (bd_isolate_summary.tsv):
  mean_depth, median_depth     over 10 kb windows
  snps, snps_per_mb            biallelic SNPs against the reference with DP >= MIN_DP and DP within
                               0.5x..2x of the median depth
  het_snps, het_fraction       SNPs whose minor-allele fraction (MAF) is 0.2 to 0.5 inclusive
  hom_alt_fraction             SNPs with alt fraction >= 0.8
  het_per_mb                   het_snps per Mb of the genome (windows with depth)
  het_windows_fraction         share of 10 kb windows holding >= 3 het SNPs
  depth_flag_seqs              sequences of >= 200 kb whose median window depth is outside 0.75..1.25 of the
                               genome median (possible copy-number change; not a ploidy call)
Per sequence (bd_depth_by_seq.tsv): median depth ratio and het SNPs per Mb.
Per window (bd_windows.tsv.gz): run, seq, start, depth_ratio, snps, het.

Heterozygosity here is allele balance in the reads, not a diploid genotype call. A reference-divergent
haploid isolate shows many hom_alt SNPs and few het SNPs; a hybrid or diploid shows het SNPs with MAF near 0.5.
Read mapping errors in repeats and paralogs also create het-like SNPs, so depth-filtered SNPs are used and
the windows table lets a reader remove outlier windows.
"""
import argparse
import collections
import gzip
import os
import statistics

MIN_DP = 10


def read_windows(path):
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, s, e, d = line.rstrip("\n").split("\t")[:4]
            rows.append((c, int(s), int(e), float(d)))
    return rows


def parse_snps(path, median_depth):
    lo, hi = 0.5 * median_depth, 2.0 * median_depth
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or "," in f[3]:
                continue  # multiallelic
            try:
                dp = int(f[5]); ad = [int(x) for x in f[6].split(",")]
            except ValueError:
                continue
            if len(ad) != 2 or dp < MIN_DP or not (lo <= dp <= hi) or sum(ad) == 0:
                continue
            yield f[0], int(f[1]), ad[1] / sum(ad)


def summarise_run(run, base, ref):
    win = read_windows(f"{base}/per_run_{ref}/{run}.regions.bed.gz")
    depths = [w[3] for w in win if w[3] > 0]
    if not depths:
        return None
    med = statistics.median(depths)
    nwin = collections.defaultdict(lambda: [0, 0])  # (chrom,start)-> [snps, het]
    per_seq = collections.defaultdict(list)
    for c, s, e, d in win:
        per_seq[c].append(d)
    snps = het = hom = 0
    seq_het = collections.Counter()
    for c, pos, af in parse_snps(f"{base}/per_run_{ref}/{run}.snps.tsv.gz", med):
        w = (c, (pos - 1) // 10000 * 10000)
        nwin[w][0] += 1
        snps += 1
        if 0.2 <= af <= 0.8:
            nwin[w][1] += 1; het += 1; seq_het[c] += 1
        elif af > 0.8:
            hom += 1
    gsize = sum(e - s for _, s, e, d in win if d > 0)
    flags, seq_rows = 0, []
    seq_len = collections.Counter()
    for c, s, e, d in win:
        seq_len[c] += e - s
    for c, ds in per_seq.items():
        r = statistics.median(ds) / med if med else 0
        L = seq_len[c]
        if L >= 200000 and not (0.75 <= r <= 1.25):
            flags += 1
        seq_rows.append((run, c, L, round(r, 3), round(seq_het[c] / (L / 1e6), 1)))
    hw = sum(1 for v in nwin.values() if v[1] >= 3)
    row = dict(run=run, mean_depth=round(statistics.mean(depths), 1), median_depth=round(med, 1), snps=snps,
               snps_per_mb=round(snps / (gsize / 1e6), 1) if gsize else 0, het_snps=het,
               het_fraction=round(het / snps, 4) if snps else 0, hom_alt_fraction=round(hom / snps, 4) if snps else 0,
               het_per_mb=round(het / (gsize / 1e6), 1) if gsize else 0,
               het_windows_fraction=round(hw / max(1, len(win)), 4), depth_flag_seqs=flags)
    wrows = [(run, c, s, round(d / med, 3), nwin[(c, s)][0], nwin[(c, s)][1]) for c, s, e, d in win]
    return row, seq_rows, wrows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True)
    ap.add_argument("--ref", default="JEL423", help="reference name used in per_run_<ref>")
    a = ap.parse_args()
    meta = {}
    for l in list(open(f"{a.base}/bd_runs.tsv"))[1:]:
        f = l.rstrip("\n").split("\t")
        meta[f[1]] = dict(group=f[0], isolate=f[3], country=f[4])
    out_sum, out_seq = [], []
    with gzip.open(f"{a.base}/bd_windows_{a.ref}.tsv.gz", "wt") as fw:
        fw.write("run\tseq\tstart\tdepth_ratio\tsnps\thet\n")
        for run, m in meta.items():
            if not (os.path.exists(f"{a.base}/per_run_{a.ref}/{run}.snps.tsv.gz") and
                    os.path.exists(f"{a.base}/per_run_{a.ref}/{run}.regions.bed.gz")):
                continue
            r = summarise_run(run, a.base, a.ref)
            if r is None:
                continue
            row, seqs, wrows = r
            row.update(m)
            out_sum.append(row); out_seq += seqs
            for w in wrows:
                fw.write("\t".join(map(str, w)) + "\n")
    cols = ["run", "group", "isolate", "country", "mean_depth", "median_depth", "snps", "snps_per_mb", "het_snps",
            "het_fraction", "hom_alt_fraction", "het_per_mb", "het_windows_fraction", "depth_flag_seqs"]
    with open(f"{a.base}/bd_isolate_summary_{a.ref}.tsv", "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in out_sum:
            fo.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(f"{a.base}/bd_depth_by_seq_{a.ref}.tsv", "w") as fo:
        fo.write("run\tseq\tlength\tdepth_ratio\thet_per_mb\n")
        for r in out_seq:
            fo.write("\t".join(map(str, r)) + "\n")
    print(len(out_sum), "isolates summarised")


if __name__ == "__main__":
    main()
