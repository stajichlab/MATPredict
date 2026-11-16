#!/usr/bin/env python3
"""Hybrid versus aneuploid pattern for chosen isolates: per-chromosome depth, allele balance and parent-of-origin window classes.

Read-only analysis. For each --run:
  JEL423 reference : per sequence (>= --min-len): depth ratio (bd_depth_by_seq_JEL423.tsv), het SNPs per Mb and the share of het sites near 1/3 and
                     near 1/2 (minor allele fraction), from per_run_JEL423/<run>.snps.tsv.gz
  CLFT044 reference: per sequence (>= --min-len): share of 50 kb windows of the parent-of-origin analysis in each class (CLFT, HET, JELHOM)
                     (parent_origin/parent_origin_windows.tsv.gz)
A first-generation hybrid should show HET windows on every sequence and an even allele balance; an aneuploid isolate shows local departures of
depth and of the 1/3 share. Isolates of one lineage should show CLFT or JELHOM windows throughout.
Output (--out): hetrich_by_sequence.tsv
"""
import argparse
import collections
import csv
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bd_aneuploidy_af import read_snps, shares  # noqa: E402


def window_classes(path, run, lengths, min_len):
    """{seq: Counter(class)} for one run from parent_origin_windows.tsv.gz, sequences of >= min_len bp."""
    out = collections.defaultdict(collections.Counter)
    with gzip.open(path, "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == run and lengths.get(f[1], 0) >= min_len:
                out[f[1]][f[7]] += 1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True)
    ap.add_argument("--run", action="append", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-len", type=int, default=500000)
    a = ap.parse_args()
    summ = {r["run"]: r for r in csv.DictReader(open(f"{a.base}/bd_isolate_summary_JEL423.tsv"), delimiter="\t")}
    dj = list(csv.DictReader(open(f"{a.base}/bd_depth_by_seq_JEL423.tsv"), delimiter="\t"))
    dc = list(csv.DictReader(open(f"{a.base}/bd_depth_by_seq_CLFT044.tsv"), delimiter="\t"))
    lj = {r["seq"]: int(r["length"]) for r in dj}
    lc = {r["seq"]: int(r["length"]) for r in dc}
    with open(a.out, "w") as fo:
        fo.write("run\treference\tseq\tlength\tdepth_ratio\the_per_mb\tshare_third\tshare_half\twindows\tshare_CLFT\tshare_HET\tshare_JELHOM\n")
        for run in a.run:
            d = float(summ[run]["median_depth"])
            snps = read_snps(f"{a.base}/per_run_JEL423/{run}.snps.tsv.gz", 0.5 * d, 2 * d)
            for r in dj:
                if r["run"] == run and lj[r["seq"]] >= a.min_len:
                    n, t, h = shares(snps.get(r["seq"], []))
                    fo.write(f"{run}\tJEL423\t{r['seq']}\t{r['length']}\t{r['depth_ratio']}\t{n / (int(r['length']) / 1e6):.1f}\t{t:.3f}\t{h:.3f}\t\t\t\t\n")
            wc = window_classes(f"{a.base}/parent_origin/parent_origin_windows.tsv.gz", run, lc, a.min_len)
            for r in dc:
                if r["run"] == run and r["seq"] in wc:
                    c = wc[r["seq"]]
                    tot = sum(c.values())
                    fo.write(f"{run}\tCLFT044\t{r['seq']}\t{r['length']}\t{r['depth_ratio']}\t\t\t\t{tot}\t{c['CLFT'] / tot:.3f}\t{c['HET'] / tot:.3f}\t{c['JELHOM'] / tot:.3f}\n")
    print(open(a.out).read())


if __name__ == "__main__":
    main()
