#!/usr/bin/env python3
"""Table of the depth-controlled site counts per isolate (output of bd_site_genotype.slurm).

Reads site_counts/<run>.counts.tsv.gz (contig, pos, n_clft, n_jel) and the isolate classes, and writes for each run: the number of divergent
sites with at least --min-depth reads, the share with a JEL fraction f <= 0.2 (clft), 0.2 < f < 0.8 (het) and >= 0.8 (jel), and the same
shares on each of the sequences of --min-len bp or more when --by-sequence is given.
Output: TSV on stdout or --out.
"""
import argparse
import collections
import csv
import gzip
import os


def shares_from_counts(path, min_depth=5):
    """-> (n_sites, {'clft': share, 'het': share, 'jel': share}, per_seq {seq: (n, shares)})."""
    tally = collections.Counter()
    per = collections.defaultdict(collections.Counter)
    with gzip.open(path, "rt") as fh:
        next(fh)
        for line in fh:
            c, p, nr, nj = line.rstrip("\n").split("\t")
            nr, nj = int(nr), int(nj)
            if nr + nj < min_depth:
                continue
            f = nj / (nr + nj)
            k = "clft" if f <= 0.2 else "jel" if f >= 0.8 else "het"
            tally[k] += 1
            per[c][k] += 1
    n = sum(tally.values())
    sh = {k: (tally[k] / n if n else float("nan")) for k in ("clft", "het", "jel")}
    ps = {}
    for c, t in per.items():
        m = sum(t.values())
        ps[c] = (m, {k: t[k] / m for k in ("clft", "het", "jel")})
    return n, sh, ps


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="site_counts folder")
    ap.add_argument("--classes", required=True)
    ap.add_argument("--summary", required=True, help="bd_isolate_summary_CLFT044.tsv (median depth)")
    ap.add_argument("--out")
    ap.add_argument("--min-depth", type=int, default=5)
    ap.add_argument("--lengths", help="bd_depth_by_seq_CLFT044.tsv; with --by-sequence-out writes per-sequence shares for sequences of >= --min-len bp")
    ap.add_argument("--by-sequence-out")
    ap.add_argument("--min-len", type=int, default=500000)
    a = ap.parse_args()
    cls = {r["run"]: r for r in csv.DictReader(open(a.classes), delimiter="\t")}
    dep = {r["run"]: r["median_depth"] for r in csv.DictReader(open(a.summary), delimiter="\t")}
    per_run = {}
    lines = ["run\tisolate\tclass\tmedian_depth_full\tsites\tshare_clft\tshare_het\tshare_jel"]
    for f in sorted(os.listdir(a.dir)):
        if not f.endswith(".counts.tsv.gz"):
            continue
        run = f[:-len(".counts.tsv.gz")]
        n, sh, ps = shares_from_counts(os.path.join(a.dir, f), a.min_depth)
        per_run[run] = ps
        lines.append(f"{run}\t{cls[run]['isolate']}\t{cls[run]['class']}\t{dep[run]}\t{n}\t{sh['clft']:.3f}\t{sh['het']:.3f}\t{sh['jel']:.3f}")
    if a.by_sequence_out and a.lengths:
        L = {}
        for r in csv.DictReader(open(a.lengths), delimiter="\t"):
            L[r["seq"]] = int(r["length"])
        with open(a.by_sequence_out, "w") as fo:
            fo.write("run\tisolate\tseq\tlength\tsites\tshare_clft\tshare_het\tshare_jel\n")
            for run, ps in per_run.items():
                for seq, (m, sh) in sorted(ps.items()):
                    if L.get(seq, 0) >= a.min_len:
                        fo.write(f"{run}\t{cls[run]['isolate']}\t{seq}\t{L[seq]}\t{m}\t{sh['clft']:.3f}\t{sh['het']:.3f}\t{sh['jel']:.3f}\n")
    text = "\n".join(lines) + "\n"
    print(text)
    if a.out:
        open(a.out, "w").write(text)


if __name__ == "__main__":
    main()
