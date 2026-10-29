#!/usr/bin/env python3
"""Group the runs of one species into clones (runs that are the same genotype) from the read sketches written by genotype_mat_one_accession.sh with SKETCH=1.

For every pair of runs: the containment of the smaller sketch in the larger (hashes seen at least MIN_ABUND times, which drops most sequencing errors), converted to an identity
estimate ANI = containment^(1/21). Two runs are in the same clone when the estimate is >= --ani (default 0.998). Clusters are made by average linkage (default; --linkage single is connected components, which chains: pilot 2026-10-27); the least
similar pair inside each is reported as `min_ani`. This reports what the reads say; it does not know the biology: sibling progeny of one cross
fall below the threshold and a clonal lineage with a few thousand SNPs can fall above a looser one. `--pairs` writes every pair identity so the distribution can be inspected; `--from-pairs` clusters it again at another threshold or linkage.

    clone_clusters.py --sigs DIR_WITH_sig.zip [--runs results_all_final.tsv --group SLUG] [--ani 0.998] [--linkage average|complete|single] [--from-pairs pairs.tsv] --out clusters.tsv [--pairs pairs.tsv]
Output columns (clusters.tsv): run, group, clone (c1, c2, ...; singletons get their own), clone_size, representative (best-covered run of the clone), min_ani, call (from --runs, if given).
`counts()` gives, per group, the MAT calls counted once per clone (a clone with mixed calls is counted as `conflict`).
"""
import argparse, csv, glob, itertools, os, sys
from collections import Counter, defaultdict

K = 21


def load(sigdir, runs):
    import sourmash
    out = {}
    for p in sorted(glob.glob(os.path.join(sigdir, "*.sig.zip"))):
        run = os.path.basename(p)[:-8]
        if runs is not None and run not in runs: continue
        for ss in sourmash.load_file_as_signatures(p, ksize=K):
            mh = ss.minhash
            out[run] = {h: a for h, a in mh.hashes.items()}
            break
    return out


def pairs(sk, min_abund):
    sets = {r: {h for h, a in d.items() if a >= min_abund} for r, d in sk.items()}
    res = []
    for a, b in itertools.combinations(sorted(sets), 2):
        A, B = sets[a], sets[b]
        n = min(len(A), len(B))
        if n < 1000: res.append((a, b, None)); continue          # too few hashes to say (low coverage or a tiny subsample)
        c = len(A & B) / n
        res.append((a, b, c ** (1.0 / K) if c > 0 else 0.0))
    return res, {r: len(s) for r, s in sets.items()}


def cluster(runs, prs, thr, linkage="average"):
    """Clusters of runs. average (default): merge the two clusters with the highest mean pair identity while it is >= thr (stops chaining); complete: lowest pair identity; single: connected components."""
    if linkage != "single":
        P = {frozenset((a, b)): v for a, b, v in prs if v is not None}
        cl = [[r] for r in runs]
        sim = lambda A, B: (lambda v: min(v) if linkage == "complete" else sum(v) / len(v))([P.get(frozenset((a, b)), 0.0) for a in A for b in B])
        while True:
            best = None
            for i in range(len(cl)):
                for j in range(i + 1, len(cl)):
                    v = sim(cl[i], cl[j])
                    if v >= thr and (best is None or v > best[0]): best = (v, i, j)
            if not best: return cl
            _, i, j = best; cl[i] += cl[j]; del cl[j]
    par = {r: r for r in runs}
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for a, b, v in prs:
        if v is not None and v >= thr: par[f(a)] = f(b)
    comp = defaultdict(list)
    for r in runs: comp[f(r)].append(r)
    return list(comp.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sigs"); ap.add_argument("--runs"); ap.add_argument("--group")
    ap.add_argument("--ani", type=float, default=0.998); ap.add_argument("--linkage", choices=["average", "complete", "single"], default="average"); ap.add_argument("--from-pairs", help="cluster a saved pairs.tsv again (no sketches needed; --sigs is ignored)"); ap.add_argument("--min-abund", type=int, default=2)
    ap.add_argument("--out", required=True); ap.add_argument("--pairs")
    a = ap.parse_args()
    calls, group = {}, a.group or ""
    sel = None
    if a.runs:
        rows = list(csv.DictReader(open(a.runs, newline="").read().replace("\r", "").splitlines(), delimiter="\t"))
        if a.group: rows = [r for r in rows if r.get("group") == a.group]
        calls = {r["run"]: r.get("call", "") for r in rows}; sel = set(calls)
    if a.from_pairs:
        prs = [(r["run_a"], r["run_b"], float(r["ani"]) if r["ani"] else None) for r in csv.DictReader(open(a.from_pairs, newline=""), delimiter="\t")]
        sk = {x: None for p in prs for x in p[:2]}; size = {x: 0 for x in sk}
    else:
        if not a.sigs: sys.exit("--sigs or --from-pairs is required")
        sk = load(a.sigs, sel)
        if len(sk) < 2: sys.exit("fewer than two sketches")
        prs, size = pairs(sk, a.min_abund)
    if a.pairs:
        with open(a.pairs, "w") as fh:
            fh.write("run_a\trun_b\tani\n"); [fh.write(f"{x}\t{y}\t{'' if v is None else f'{v:.5f}'}\n") for x, y, v in prs]
    pv = {frozenset((x, y)): v for x, y, v in prs}
    comps = sorted(cluster(sorted(sk), prs, a.ani, a.linkage), key=lambda c: (-len(c), c[0]))
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["run", "group", "clone", "clone_size", "representative", "min_ani", "call"])
        for i, c in enumerate(comps, 1):
            rep = max(c, key=lambda r: size[r]); vals = [pv[frozenset(p)] for p in itertools.combinations(c, 2) if pv.get(frozenset(p)) is not None]
            for r in c: w.writerow([r, group, f"c{i}", len(c), rep, f"{min(vals):.4f}" if vals else "", calls.get(r, "")])
    sizes = Counter(len(c) for c in comps)
    print(f"{len(sk)} runs -> {len(comps)} clones; sizes {dict(sorted(sizes.items()))}")


def counts(clusters_tsv):
    """{group: Counter of calls per clone} from a clusters.tsv written with --runs."""
    by = defaultdict(lambda: defaultdict(set))
    for r in csv.DictReader(open(clusters_tsv, newline=""), delimiter="\t"):
        if r["call"] in ("MAT1-1", "MAT1-2", "both"): by[r["group"]][r["clone"]].add(r["call"])
    return {g: Counter(next(iter(s)) if len(s) == 1 else "conflict" for s in cl.values()) for g, cl in by.items()}


if __name__ == "__main__":
    main()
