#!/usr/bin/env python3
"""Depth of a read subsample along each idiomorph of a panel, in windows: is a `both` or `none` call real?

    depth_profile.py --reads sub.fastq --panel-dir DIR --prefix PFX --out profile.tsv [--window 500] [--idiomorphs "MAT1-1 MAT1-2"]
For every window of every panel sequence: mean count of its k-mers in the reads, and the mean over the idiomorph-specific k-mers only. A real idiomorph
is covered along its whole length at a depth near the single-copy flank; a trace (cross-talk, contamination) covers short stretches or a few k-mers.
Panels with several records per idiomorph are profiled record by record.
"""
import argparse, collections, sys
from MATPredict.reads.panel import Panel
from MATPredict.reads.kmers import iter_fastq_reads, revcomp


def records(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name: yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        elif line.strip(): seq.append(line.strip())
    if name: yield name, "".join(seq)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--reads", required=True); ap.add_argument("--panel-dir", required=True); ap.add_argument("--prefix", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--window", type=int, default=500); ap.add_argument("--idiomorphs", default="MAT1-1 MAT1-2")
    a = ap.parse_args()
    ids = a.idiomorphs.split(); k = 31
    fastas = {i: f"{a.panel_dir}/{a.prefix}_{i}.fasta" for i in ids}
    p = Panel.from_fastas(fastas, k=k, min_run=100)
    table, classes = p.lookup()
    cnt = collections.Counter(); n = 0
    for r in iter_fastq_reads(a.reads):
        n += 1
        for i in range(len(r) - k + 1):
            km = r[i:i + k]
            if km in table: cnt[km] += 1
    with open(a.out, "w") as fh:
        fh.write("idiomorph\trecord\twindow_start\tmean_depth_all\tmean_depth_unique\tn_unique\tn_shared\n")
        for idio, path in fastas.items():
            for rec, seq in records(path):
                for st in range(0, max(1, len(seq) - k + 1), a.window):
                    tot = uni = nu = ns = 0; m = 0
                    for j in range(st, min(st + a.window, len(seq) - k + 1)):
                        km = seq[j:j + k]; rc = revcomp(km)
                        c = cnt.get(km, 0) + (cnt.get(rc, 0) if rc != km else 0)
                        m += 1; tot += c
                        hit = table.get(km) or table.get(rc)
                        if hit and hit[0] == idio: uni += c; nu += 1
                        elif hit and hit[0] == "shared": ns += 1
                    fh.write(f"{idio}\t{rec}\t{st}\t{tot / max(m, 1):.2f}\t{uni / max(nu, 1):.2f}\t{nu}\t{ns}\n")
    print("reads", n, file=sys.stderr)


if __name__ == "__main__":
    main()
