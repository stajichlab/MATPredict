#!/usr/bin/env python3
"""Compare self-map depth and read mismatch of an isolate with the JEL-base share of its own assembly, per window.

Read-only analysis. Inputs: bd_self_map_windows.slurm output (depth_windows.tsv, mismatch_windows.tsv; windows on the isolate's own
assembly), the PAF of that assembly (query) against CLFT044 (target), and site_states.tsv.gz from bd_ref_lineage_sites.py (CLFT044
coordinates; one state column per assembly). Each site is moved to assembly coordinates by linear interpolation inside its primary
PAF block (an approximation; error is far smaller than a 100 kb window). Per window: n sites with state jel or clft in the assembly,
JEL share = jel / (jel + clft), median depth, edit distance per aligned base.
Usage: bd_selfmap_vs_jel_share.py --paf P.paf.gz --states site_states.tsv.gz --column CLFT067 --depth D.tsv --mismatch M.tsv --out OUT.tsv [--min-block 5000]
"""
import argparse
import bisect
import collections
import gzip
import sys


def open_text(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def load_blocks(paf, min_block):
    """-> {target contig: sorted [(ts, te, qname, qs, qe, strand)]} for primary blocks >= min_block (0-based half-open)."""
    out = collections.defaultdict(list)
    with open_text(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or int(f[10]) < min_block:
                continue
            if any(t == "tp:A:S" for t in f[12:]):
                continue
            out[f[5]].append((int(f[7]), int(f[8]), f[0], int(f[2]), int(f[3]), f[4]))
    for v in out.values():
        v.sort()
    return out


def lift(blocks, contig, pos0):
    """Target 0-based position -> (query name, query position) or None when no block, or more than one block, covers it."""
    v = blocks.get(contig)
    if not v:
        return None
    i = bisect.bisect_right(v, (pos0, float("inf"))) 
    hits = [b for b in v[max(0, i - 50):i] if b[0] <= pos0 < b[1]]
    if len(hits) != 1:
        return None
    ts, te, q, qs, qe, st = hits[0]
    frac = (pos0 - ts) / max(1, te - ts)
    qp = qs + frac * (qe - qs) if st == "+" else qe - frac * (qe - qs)
    return q, int(qp)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--states", required=True)
    ap.add_argument("--column", required=True)
    ap.add_argument("--depth", required=True)
    ap.add_argument("--mismatch", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=int, default=100000)
    ap.add_argument("--min-block", type=int, default=5000)
    a = ap.parse_args()
    blocks = load_blocks(a.paf, a.min_block)
    cnt = collections.defaultdict(lambda: [0, 0])
    with open_text(a.states) as fh:
        head = next(fh).rstrip("\n").split("\t")
        ci = head.index(a.column)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            st = f[ci]
            if st not in ("jel", "clft"):
                continue
            r = lift(blocks, f[0], int(f[1]) - 1)
            if r is None:
                continue
            cnt[(r[0], r[1] // a.window)][st == "jel"] += 1
    depth = {}
    with open(a.depth) as fh:
        next(fh)
        for line in fh:
            c, s, e, d = line.rstrip("\n").split("\t")
            depth[(c, int(s) // a.window)] = float(d)
    mism = {}
    with open(a.mismatch) as fh:
        next(fh)
        for line in fh:
            c, s, n, bp, ed = line.rstrip("\n").split("\t")
            mism[(c, int(s) // a.window)] = (int(bp), int(ed))
    with open(a.out, "w") as fo:
        fo.write("contig\twindow_start\tn_sites\tjel_share\tdepth\tedit_per_base\n")
        for k in sorted(depth):
            c = cnt.get(k, [0, 0])
            n = c[0] + c[1]
            bp, ed = mism.get(k, (0, 0))
            fo.write(f"{k[0]}\t{k[1] * a.window}\t{n}\t{(c[1] / n if n else float('nan')):.4f}\t{depth[k]:.2f}\t{(ed / bp if bp else float('nan')):.5f}\n")


if __name__ == "__main__":
    sys.exit(main())
