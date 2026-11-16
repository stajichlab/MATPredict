#!/usr/bin/env python3
"""Allele balance on flagged chromosomes: does a depth gain look like an extra chromosome copy?

For each flag of bd_aneuploidy.py (aneuploidy_flags.tsv) with a normalised ratio >= --gain or <= --loss, take the heterozygous SNPs of
that isolate (bcftools table per_run_<ref>/<run>.snps.tsv.gz; biallelic, DP >= 10 and within 0.5x to 2x of the isolate median depth, as in
bd_summarize.py) on the flagged sequence and on all other sequences of 500 kb or more, fold the alternative-allele fraction to the minor
allele (<= 0.5), and report the share of het sites (minor fraction 0.15 to 0.5) in three bins: near 1/3 (0.27 to 0.40), near 1/2 (0.43 to 0.50)
and other. A trisomic chromosome of a diploid should move sites from 1/2 to 1/3; a monosomic one removes het sites (het per Mb falls).
With --genomewide the per-isolate shares over all sequences of at least --min-len and a baseline label (third-dominated if the 1/3 share exceeds the
1/2 share by 0.1, half-dominated in the opposite case, mixed otherwise) are written too: a first look at the baseline ploidy state of each isolate.
Output: aneuploidy_af.tsv (run, isolate, seq, normalised, het_sites_flagged, share_third_flagged, share_half_flagged, het_sites_rest,
share_third_rest, share_half_rest, het_per_mb_flagged, het_per_mb_rest).
Limit: low-depth isolates give wide allele-fraction noise (a 20x site has binomial sd of about 0.11), so a shift from 0.5 to 0.33 needs many sites.
"""
import argparse
import csv
import gzip
import statistics
import collections


def fold(af):
    return min(af, 1 - af)


def shares(mafs):
    """(n, share near 1/3, share near 1/2) among minor-allele fractions of het sites (0.15 to 0.5)."""
    het = [m for m in mafs if 0.15 <= m <= 0.5]
    if not het:
        return 0, float("nan"), float("nan")
    third = sum(0.27 <= m <= 0.40 for m in het) / len(het)
    half = sum(0.43 <= m <= 0.50 for m in het) / len(het)
    return len(het), third, half


def read_snps(path, lo, hi):
    out = collections.defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or "," in f[3]:
                continue
            try:
                dp = int(f[5]); ad = [int(x) for x in f[6].split(",")]
            except ValueError:
                continue
            if len(ad) == 2 and sum(ad) and dp >= 10 and lo <= dp <= hi:
                out[f[0]].append(fold(ad[1] / sum(ad)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flags", required=True)
    ap.add_argument("--classes", help="isolate_classes.tsv (only for --genomewide labels)")
    ap.add_argument("--base", required=True)
    ap.add_argument("--ref", default="JEL423")
    ap.add_argument("--lengths", required=True, help="bd_depth_by_seq_<ref>.tsv (for sequence lengths)")
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gain", type=float, default=1.3)
    ap.add_argument("--loss", type=float, default=0.7)
    ap.add_argument("--min-len", type=int, default=500000)
    ap.add_argument("--genomewide", help="also write per-isolate genome-wide allele-balance shares (all isolates with a snp table) to this TSV")
    a = ap.parse_args()
    L = {r["seq"]: int(r["length"]) for r in csv.DictReader(open(a.lengths), delimiter="\t")}
    depth = {r["run"]: float(r["median_depth"]) for r in csv.DictReader(open(a.summary), delimiter="\t")}
    big = {s for s, l in L.items() if l >= a.min_len}
    rows = []
    for f in csv.DictReader(open(a.flags), delimiter="\t"):
        nv = float(f["normalised"])
        if not (nv >= a.gain or nv <= a.loss):
            continue
        d = depth[f["run"]]
        snps = read_snps(f"{a.base}/per_run_{a.ref}/{f['run']}.snps.tsv.gz", 0.5 * d, 2 * d)
        fl = snps.get(f["seq"], [])
        rest = [m for s, v in snps.items() if s in big and s != f["seq"] for m in v]
        n1, t1, h1 = shares(fl)
        n2, t2, h2 = shares(rest)
        lenf = L[f["seq"]] / 1e6
        lenr = sum(L[s] for s in big if s != f["seq"]) / 1e6
        rows.append((f["run"], f["isolate"], f["seq"], f["normalised"], n1, round(t1, 3), round(h1, 3), n2, round(t2, 3), round(h2, 3),
                     round(n1 / lenf, 1), round(n2 / lenr, 1)))
    if a.genomewide:
        import os
        cl = {r["run"]: r for r in csv.DictReader(open(a.classes), delimiter="\t")} if a.classes else {}
        with open(a.genomewide, "w") as fg:
            fg.write("run\tisolate\tclass\tmedian_depth\the_sites\tshare_third\tshare_half\tbaseline\n")
            for run, d in depth.items():
                path = f"{a.base}/per_run_{a.ref}/{run}.snps.tsv.gz"
                if not os.path.exists(path):
                    continue
                snps = read_snps(path, 0.5 * d, 2 * d)
                allm = [m for sq, v in snps.items() if sq in big for m in v]
                n, t, h = shares(allm)
                base = "NA" if n < 500 else ("third-dominated" if t > h + 0.1 else "half-dominated" if h > t + 0.1 else "mixed")
                fg.write(f"{run}\t{cl.get(run, {}).get('isolate', '')}\t{cl.get(run, {}).get('class', '')}\t{d}\t{n}\t{round(t, 3)}\t{round(h, 3)}\t{base}\n")
    with open(a.out, "w") as fo:
        fo.write("run\tisolate\tseq\tnormalised\the_flagged\tshare_third_flagged\tshare_half_flagged\the_rest\tshare_third_rest\tshare_half_rest\the_per_mb_flagged\the_per_mb_rest\n")
        for r in sorted(rows, key=lambda r: -abs(float(r[3]) - 1)):
            fo.write("\t".join(map(str, r)) + "\n")
    print(open(a.out).read())


if __name__ == "__main__":
    main()
