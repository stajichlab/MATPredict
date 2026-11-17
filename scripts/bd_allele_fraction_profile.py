#!/usr/bin/env python3
"""JEL423-allele fraction at lineage-diagnostic sites, per isolate: histogram and genome windows.

Read-only analysis. Input is the per-run read counts at the divergent sites (bd_site_allele_counts.py output: contig, pos, n_clft, n_jel)
and a list of sites to keep (for example diagnostic_sites_no067.tsv: contig, pos, JEL423 base).
For each isolate, sites with at least --min-depth reads are used. f = n_jel / (n_clft + n_jel).
  hist     counts of sites in ten bins of f (0-0.1 ... 0.9-1.0) and the share of sites with 0.3 <= f <= 0.7 ("mid")
  windows  per window of --window bp on a contig (contigs shorter than the window are one window): sites, median f, mid share.
A heterozygous isolate gives a peak near f = 0.5 (a binomial spread that depends on depth); a haploid or homozygous one gives peaks near 0 and 1.
At low depth the binomial spread is wide: with depth d, the share of true f = 0.5 sites inside 0.3-0.7 is reported as 'expected_mid_if_het'
(computed from the depth of each site; a model value, not data).
Usage: bd_allele_fraction_profile.py --sites SITES.tsv --counts NAME=PATH [...] --out DIR [--min-depth 8] [--window 100000]
"""
import argparse
import collections
import gzip
import math
import os
import statistics
import sys


def open_text(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def binom_mid(n, lo=0.3, hi=0.7, p=0.5):
    """P(lo <= k/n <= hi) for k ~ Binomial(n, p)."""
    return sum(math.comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(n + 1) if lo <= k / n <= hi)


def load_sites(path):
    return {(l.split("\t")[0], int(l.split("\t")[1])) for l in open_text(path) if l.strip()}


def profile(counts_path, sites, min_depth, window):
    """-> (hist list of 10, n, mid_n, expected_mid_sum, windows {(contig, idx): [f, ...]})."""
    hist = [0] * 10
    n = mid = 0
    exp = 0.0
    cache = {}
    win = collections.defaultdict(list)
    with open_text(counts_path) as fh:
        next(fh)
        for line in fh:
            c, p, a, b = line.rstrip("\n").split("\t")
            key = (c, int(p))
            if key not in sites:
                continue
            a, b = int(a), int(b)
            d = a + b
            if d < min_depth:
                continue
            f = b / d
            hist[min(9, int(f * 10))] += 1
            n += 1
            mid += 0.3 <= f <= 0.7
            if d not in cache:
                cache[d] = binom_mid(d)
            exp += cache[d]
            win[(c, key[1] // window)].append(f)
    return hist, n, mid, exp, win


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--counts", action="append", required=True, help="NAME=PATH")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-depth", type=int, default=8)
    ap.add_argument("--window", type=int, default=100000)
    ap.add_argument("--min-window-sites", type=int, default=20)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sites = load_sites(a.sites)
    with open(os.path.join(a.out, "allele_fraction_hist.tsv"), "w") as fh, open(os.path.join(a.out, "allele_fraction_windows.tsv"), "w") as fw:
        fh.write("isolate\tn_sites\t" + "\t".join(f"f{i / 10:.1f}" for i in range(10)) + "\tmid_share\texpected_mid_if_het\n")
        fw.write("isolate\tcontig\twindow_start\tn_sites\tmedian_f\tmid_share\n")
        for spec in a.counts:
            name, path = spec.split("=", 1)
            hist, n, mid, exp, win = profile(path, sites, a.min_depth, a.window)
            ms = mid / n if n else float("nan")
            es = exp / n if n else float("nan")
            fh.write(f"{name}\t{n}\t" + "\t".join(map(str, hist)) + f"\t{ms:.3f}\t{es:.3f}\n")
            for (c, i), v in sorted(win.items()):
                if len(v) >= a.min_window_sites:
                    fw.write(f"{name}\t{c}\t{i * a.window}\t{len(v)}\t{statistics.median(v):.3f}\t{sum(0.3 <= x <= 0.7 for x in v) / len(v):.3f}\n")
    print(open(os.path.join(a.out, "allele_fraction_hist.tsv")).read())


if __name__ == "__main__":
    sys.exit(main())
