#!/usr/bin/env python3
"""Isolate-level SNP ordination (PCA) of the Bd isolates from the SNP tables against one reference.

Input  --base       results folder with per_run_<ref>/<run>.snps.tsv.gz (bcftools table: CHROM POS REF ALT QUAL DP AD) and
                    bd_isolate_summary_<ref>.tsv (median_depth); --classes isolate_classes.tsv for labels
Method : biallelic SNPs on sequences of at least --min-len with DP >= 10 and within 0.5x to 2x of the isolate median depth. To keep memory
         small only sites with POS % --thin == 0 are used. A site enters the matrix when at least --min-called isolates carry it and at
         most --max-called share of isolates carry it (the fixed differences between the reference and almost all isolates carry no
         structure). Entry = alternative allele fraction in that isolate; an isolate without a SNP row at the site gets 0 (the reference allele;
         a missing row may also mean no coverage, so the matrix is biased toward zero at low depth). Columns are mean-centred and the
         isolate x isolate PCA comes from the SVD.
Output : --out pca_scores.tsv (run, isolate, class, PC1..PC4, sites_called), pca_variance.tsv.
Limit  : alt fraction is not a genotype call; low-depth isolates are noisy; the PCA shows relatedness in allele-fraction space, not a tree.
"""
import argparse
import collections
import csv
import gzip
import os


def thinned_sites(path, lo, hi, thin, bigset):
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or "," in f[3] or f[0] not in bigset:
                continue
            pos = int(f[1])
            if pos % thin:
                continue
            try:
                dp = int(f[5]); ad = [int(x) for x in f[6].split(",")]
            except ValueError:
                continue
            if len(ad) == 2 and sum(ad) and dp >= 10 and lo <= dp <= hi:
                out[(f[0], pos)] = ad[1] / sum(ad)
    return out


def main():
    import numpy as np
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True)
    ap.add_argument("--ref", default="CLFT044")
    ap.add_argument("--classes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-len", type=int, default=500000)
    ap.add_argument("--thin", type=int, default=10)
    ap.add_argument("--min-called", type=int, default=3)
    ap.add_argument("--max-called", type=float, default=0.9)
    a = ap.parse_args()
    summ = {r["run"]: float(r["median_depth"]) for r in csv.DictReader(open(f"{a.base}/bd_isolate_summary_{a.ref}.tsv"), delimiter="\t")}
    cls = {r["run"]: r for r in csv.DictReader(open(a.classes), delimiter="\t")}
    runs = [r for r in summ if os.path.exists(f"{a.base}/per_run_{a.ref}/{r}.snps.tsv.gz") and r in cls]
    lens = collections.defaultdict(int)
    with open(f"{a.base}/bd_depth_by_seq_{a.ref}.tsv") as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            lens[f[1]] = int(f[2])
    big = {s for s, l in lens.items() if l >= a.min_len}
    data = {}
    for r in runs:
        d = summ[r]
        data[r] = thinned_sites(f"{a.base}/per_run_{a.ref}/{r}.snps.tsv.gz", 0.5 * d, 2 * d, a.thin, big)
    count = collections.Counter(s for v in data.values() for s in v)
    keep = [s for s, n in count.items() if n >= a.min_called and n <= a.max_called * len(runs)]
    idx = {s: i for i, s in enumerate(keep)}
    X = np.zeros((len(runs), len(keep)))
    for i, r in enumerate(runs):
        for s, v in data[r].items():
            j = idx.get(s)
            if j is not None:
                X[i, j] = v
    X -= X.mean(axis=0)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    var = S ** 2 / (S ** 2).sum()
    with open(a.out + "/pca_scores.tsv", "w") as fo:
        fo.write("run\tisolate\tclass\tgroup\tPC1\tPC2\tPC3\tPC4\tsites_called\n")
        for i, r in enumerate(runs):
            sc = U[i, :4] * S[:4]
            fo.write(f"{r}\t{cls[r]['isolate']}\t{cls[r]['class']}\t{cls[r]['group']}\t" + "\t".join(f"{x:.2f}" for x in sc) + f"\t{len(data[r])}\n")
    with open(a.out + "/pca_variance.tsv", "w") as fo:
        fo.write("PC\tvariance_share\n")
        for k in range(6):
            fo.write(f"PC{k + 1}\t{var[k]:.4f}\n")
    print(len(runs), "isolates;", len(keep), "sites; variance PC1-4:", [round(float(x), 3) for x in var[:4]])


if __name__ == "__main__":
    main()
