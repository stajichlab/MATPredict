#!/usr/bin/env python3
"""Clone groups of every species of the reads survey from the worker's read sketches, by "solid containment" ANI (the distance of
datasets/reads/sketch_store/sketch_query.py, validated 2026-11-07 on 136 A. fumigatus strains against SNP calls).

Solid containment: a hash counts in run a's denominator only if its abundance in a is at least max(3, 0.25 x a's median abundance),
and as shared if run b has it at all. Sequencing errors (low abundance) leave the denominator, and coverage gaps in b do not count as
differences. ANI(a,b) = mean over both directions of containment^(1/21). Sketches (k=21, scaled 500, abundances; written by
genotype_mat_one_accession.sh with SKETCH=1 from the run's subsample) are downsampled to scaled 1000 first, as in sketch_query.py.

Clusters: average linkage (scipy) on 1-ANI, cut at each of THRESHOLDS. The site builder (scripts/build_reads_site.py) reads the
clone ids and applies its plateau rule; this script only records what the reads say.

    clone_solid.py --pilot PILOT_DIR --sigs DIR [--sigs DIR ...] --out OUT_DIR [--threads N] [--only SLUG ...]

PILOT_DIR: the committed wave folders (wave*/results_all_final.tsv give run -> species group). DIR: folders holding <run>.sig.zip
(searched in the order given; the first copy of a run wins). Writes per group OUT_DIR/<slug>.clones.tsv (one row per sketched run) and
OUT_DIR/same_strain_pairs.tsv (slug, strain, run_a, run_b, ani of every same-strain pair) and OUT_DIR/summary.tsv (one row per group, with a same-strain check: pairs of runs whose strain names match, their median ANI and how
many reach 0.9995 and 0.999). Columns of <slug>.clones.tsv:
  run, hashes, solid_hashes, median_abund (k-mer abundance: a coverage proxy for the subsample), nearest_run, nearest_ani,
  clone_<t> for each threshold t (c1 = the largest group), min_ani_<t> (least similar pair inside that run's group; blank for a singleton).
"""
import argparse, csv, glob, os, re, sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

K, SCALED, SOLID, MIN_ABUND = 21, 1000, 0.25, 2
MAXH = 2 ** 64
THRESHOLDS = (0.9999, 0.9995, 0.999, 0.998)
LOW_COVERAGE = 8          # median k-mer abundance below this: too few solid hashes to separate clones (sketch_query.py qc flag)


def downsample(d, scaled=SCALED, min_abund=MIN_ABUND):
    """{hash: abundance} kept to hashes below 2^64/scaled (FracMinHash downsampling) and seen at least min_abund times. The worker's sketches keep
    k-mers seen once (mostly sequencing errors, often more numerous than the genome's k-mers, so they pull the median abundance to 1-5); the sketch
    store drops them when it sketches (sketch_one.sh), and the solid threshold and the coverage proxy were validated that way."""
    lim = MAXH // scaled
    return {h: v for h, v in d.items() if h < lim and v >= min_abund}


def solid(d):
    """(median abundance, set of solid hashes) of one sketch: abundance >= max(3, 0.25 x median)."""
    if not d: return 0.0, set()
    med = float(np.median(list(d.values())))
    thr = max(3.0, SOLID * med)
    return med, {h for h, v in d.items() if v >= thr}


def ani_matrix(sketches):
    """Solid-containment ANI between every pair of sketches (list of {hash: abundance}); returns an n x n array, diagonal 1."""
    n = len(sketches)
    sol = [solid(d)[1] for d in sketches]
    alls = [set(d) for d in sketches]
    cnt = defaultdict(int)
    for s in alls:
        for h in s: cnt[h] += 1
    cols = {h: i for i, h in enumerate(h for h, c in cnt.items() if c >= 2)}
    A = np.zeros((n, len(cols)), dtype=np.float32); S = np.zeros((n, len(cols)), dtype=np.float32)
    for i in range(n):
        A[i, [cols[h] for h in alls[i] if h in cols]] = 1
        S[i, [cols[h] for h in sol[i] if h in cols]] = 1
    I = S @ A.T                                                 # |solid(a) ∩ all(b)|
    C = I / np.maximum(np.array([len(s) for s in sol], dtype=float), 1)[:, None]
    C = np.clip(C, 1e-12, 1) ** (1 / K)
    M = 0.5 * (C + C.T); np.fill_diagonal(M, 1.0)
    return M


def clusters(M, thresholds=THRESHOLDS):
    """{t: list of clone labels "c1".. (c1 the largest; ties by first member)} from average linkage on 1-M."""
    n = len(M)
    if n == 1: return {t: ["c1"] for t in thresholds}
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    D = np.clip(1 - M.astype(np.float64), 0, None); np.fill_diagonal(D, 0)
    Z = linkage(squareform(D, checks=False), method="average")
    out = {}
    for t in thresholds:
        raw = fcluster(Z, 1 - t + 1e-12, criterion="distance")
        members = defaultdict(list)
        for i, c in enumerate(raw): members[c].append(i)
        order = sorted(members.values(), key=lambda m: (-len(m), m[0]))
        lab = [""] * n
        for j, m in enumerate(order, 1):
            for i in m: lab[i] = f"c{j}"
        out[t] = lab
    return out


def group_min(M, labels):
    """Per run, the least similar pair inside its group (None for a singleton)."""
    idx = defaultdict(list)
    for i, c in enumerate(labels): idx[c].append(i)
    res = [None] * len(labels)
    for m in idx.values():
        if len(m) < 2: continue
        v = float(M[np.ix_(m, m)][np.triu_indices(len(m), 1)].min())
        for i in m: res[i] = v
    return res


NOT_A_NAME = {"", "-", "na", "n/a", "none", "missing", "not applicable", "not collected", "unknown", "not provided"}


def strain_key(name):
    """A strain name normalised for the same-strain check (case, spaces, punctuation); "" for a placeholder."""
    k = re.sub(r"[^a-z0-9]", "", (name or "").lower())
    return "" if (name or "").strip().lower() in NOT_A_NAME or len(k) < 3 else k


STRAINS = {}


def species_groups(pilot):
    """{slug: [runs]} from every committed wave table (results_all_final.tsv), runs in first-seen order, no duplicates.
    Also fills STRAINS {run: strain name} for the same-strain check."""
    g = defaultdict(list); seen = set()
    for w in sorted(glob.glob(os.path.join(pilot, "wave*"))):
        f = os.path.join(w, "results_all_final.tsv")
        if not (re.fullmatch(r"wave\d+[a-z_]*", os.path.basename(w)) and os.path.exists(f)): continue
        for r in csv.DictReader(open(f, newline=""), delimiter="\t"):
            if r.get("run") and r.get("group") and r["run"] not in seen:
                seen.add(r["run"]); g[r["group"]].append(r["run"]); STRAINS[r["run"]] = r.get("strain", "")
    return g


def load_sig(path):
    import sourmash
    for ss in sourmash.load_file_as_signatures(path, ksize=K):
        return downsample(dict(ss.minhash.hashes))
    return {}


def same_strain(M, runs, strains):
    """ANI of every pair of runs whose strain names match (after strain_key): the same strain submitted twice, a natural check of
    what 'one clone' looks like at this coverage (some are lab lines that differ, some names are reused, so it is a guide, not truth)."""
    by = defaultdict(list)
    for i, r in enumerate(runs):
        k = strain_key(strains.get(r, ""))
        if k: by[k].append(i)
    return [float(M[i, j]) for m in by.values() for x, i in enumerate(m) for j in m[x + 1:]]


def same_strain_rows(M, runs, strains):
    """(strain, run_a, run_b, ani) for every same-strain pair (the detail behind same_strain), for OUT_DIR/same_strain_pairs.tsv."""
    by = defaultdict(list)
    for i, r in enumerate(runs):
        k = strain_key(strains.get(r, ""))
        if k: by[k].append(i)
    return [(strains[runs[m[0]]], runs[i], runs[j], float(M[i, j])) for m in by.values() for x, i in enumerate(m) for j in m[x + 1:]]


def write_group(slug, runs, sketches, out, strains=None):
    M = ani_matrix(sketches)
    labs = clusters(M)
    mins = {t: group_min(M, labs[t]) for t in THRESHOLDS}
    stats = [solid(d) for d in sketches]
    n = len(runs)
    with open(os.path.join(out, slug + ".clones.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run", "hashes", "solid_hashes", "median_abund", "nearest_run", "nearest_ani"] + [f"clone_{t}" for t in THRESHOLDS] + [f"min_ani_{t}" for t in THRESHOLDS])
        for i, r in enumerate(runs):
            if n > 1:
                row = M[i].copy(); row[i] = -1; j = int(row.argmax()); near, nani = runs[j], f"{M[i, j]:.6f}"
            else:
                near, nani = "", ""
            w.writerow([r, len(sketches[i]), len(stats[i][1]), f"{stats[i][0]:.1f}", near, nani] + [labs[t][i] for t in THRESHOLDS]
                       + ["" if mins[t][i] is None else f"{mins[t][i]:.6f}" for t in THRESHOLDS])
    iu = np.triu_indices(n, 1)
    between = M[iu] if n > 1 else np.array([])
    meds = [s[0] for s in stats]
    ss = same_strain(M, runs, strains or {})
    if ss:
        with open(os.path.join(out, "same_strain_pairs.tsv"), "a", newline="") as fh:
            for st, ra, rb, v in same_strain_rows(M, runs, strains): fh.write(f"{slug}\t{st}\t{ra}\t{rb}\t{v:.6f}\n")
    return {"slug": slug, "runs": n,
            "same_strain_pairs": len(ss), "same_strain_median": f"{np.median(ss):.6f}" if ss else "",
            "same_strain_ge_0.9995": sum(v >= 0.9995 for v in ss), "same_strain_ge_0.999": sum(v >= 0.999 for v in ss),
            "median_abund": f"{np.median(meds):.1f}", "low_coverage_runs": sum(m < LOW_COVERAGE for m in meds),
            "between_median": f"{np.median(between):.6f}" if len(between) else "", "between_p99": f"{np.quantile(between, .99):.6f}" if len(between) else "",
            **{f"clones_{t}": len(set(labs[t])) for t in THRESHOLDS}}


def _load(arg):
    run, path = arg
    try: return run, load_sig(path)
    except Exception as e:
        print(f"unreadable sketch {path}: {e}", file=sys.stderr); return run, None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pilot", required=True); ap.add_argument("--sigs", action="append", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4))); ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    where = {}
    for d in a.sigs:
        for p in glob.glob(os.path.join(d, "*.sig.zip")):
            r = os.path.basename(p)[:-len(".sig.zip")]
            if "." not in r: where.setdefault(r, p)            # <run>.a.sig.zip halves are not runs
    groups = species_groups(a.pilot)
    rows = []
    with open(os.path.join(a.out, "same_strain_pairs.tsv"), "w") as fh: fh.write("slug\tstrain\trun_a\trun_b\tani\n")
    for slug in sorted(groups):
        if a.only and slug not in a.only: continue
        have = [r for r in groups[slug] if r in where]
        if not have: continue
        with Pool(a.threads) as p:
            got = [(r, d) for r, d in p.map(_load, [(r, where[r]) for r in have]) if d]
        if not got: continue
        s = write_group(slug, [r for r, _ in got], [d for _, d in got], a.out, STRAINS)
        s["typed_runs"] = len(groups[slug]); rows.append(s)
        print(f"{slug}: {s['runs']} of {len(groups[slug])} runs sketched; clones at 0.9995 {s['clones_0.9995']}; median abundance {s['median_abund']}", flush=True)
    cols = ["slug", "typed_runs", "runs", "median_abund", "low_coverage_runs", "between_median", "between_p99"] + [f"clones_{t}" for t in THRESHOLDS] + \
        ["same_strain_pairs", "same_strain_median", "same_strain_ge_0.9995", "same_strain_ge_0.999"]
    with open(os.path.join(a.out, "summary.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
