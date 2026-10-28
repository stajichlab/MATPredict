#!/usr/bin/env python3
"""Same-strain ceiling for clone_clusters.py: for every run sketched with SKETCH_SPLIT=1, the estimated identity between its two halves (two samples of one genome, each at half the coverage of the run),
beside the identities between different runs of the species. If two runs are as similar as two halves of one run, the sketches cannot tell them apart from the same strain.

    clone_calibration.py --sigs DIR --pairs pairs_<species>.tsv [--species-runs RUN,RUN,...] --out calibration.tsv
Output: one row per run with halves_ani, then a summary on stdout: the halves' range and how many between-run pairs reach the lowest halves value.
"""
import argparse, csv, glob, os, statistics
import clone_clusters as cc


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0]); ap.add_argument("--sigs", required=True); ap.add_argument("--pairs", help="pairs.tsv from clone_clusters.py (--pairs)", action="append", default=[])
    ap.add_argument("--out", required=True); ap.add_argument("--min-abund", type=int, default=2)
    a = ap.parse_args()
    runs = sorted({os.path.basename(p)[:-len(".a.sig.zip")] for p in glob.glob(os.path.join(a.sigs, "*.a.sig.zip")) if os.path.exists(p.replace(".a.sig.zip", ".b.sig.zip"))})
    rows = []
    import sourmash
    def one(p):
        for ss in sourmash.load_file_as_signatures(p, ksize=cc.K): return {h: n for h, n in ss.minhash.hashes.items()}
    for r in runs:
        A, B = one(os.path.join(a.sigs, r + ".a.sig.zip")), one(os.path.join(a.sigs, r + ".b.sig.zip"))
        prs, size = cc.pairs({r + ".a": A, r + ".b": B}, a.min_abund)
        rows.append((r, prs[0][2], size[r + ".a"], size[r + ".b"]))
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["run", "halves_ani", "hashes_a", "hashes_b"])
        for r, v, x, y in rows: w.writerow([r, "" if v is None else f"{v:.5f}", x, y])
    vals = sorted(v for _, v, _, _ in rows if v is not None)
    if not vals: print("no halves with enough hashes"); return
    print(f"{len(vals)} runs; halves ANI min {vals[0]:.4f} median {statistics.median(vals):.4f} max {vals[-1]:.4f}")
    for pf in a.pairs:
        between = [float(r["ani"]) for r in csv.DictReader(open(pf), delimiter="\t") if r["ani"]]
        print(f"{os.path.basename(pf)}: {len(between)} pairs; {sum(1 for v in between if v >= vals[0])} reach the lowest halves value {vals[0]:.4f}; {sum(1 for v in between if v >= statistics.median(vals))} reach the median {statistics.median(vals):.4f}")


if __name__ == "__main__":
    main()
