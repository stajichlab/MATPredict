#!/usr/bin/env python3
"""Cross-isolate depth scan for presence/absence regions that could hold a mating-type idiomorph (Bd).

Idea (not validated for Bd): in a diploid population with two idiomorphs, a window inside the idiomorph shows
depth ratio near 1 in isolates with two copies of the reference idiomorph, near 0.5 in hemizygous isolates (one
copy of the reference idiomorph, one copy of the other) and near 0 in isolates with none. Hemizygous windows should
carry about no het SNPs because only one copy aligns.

Input  : --base results folder with bd_windows_<ref>.tsv.gz (10 kb windows: run seq start depth_ratio snps het),
         bd_isolate_summary_<ref>.tsv (median_depth), bd_runs.tsv; --paf JEL423_vs_CLFT044.paf.gz; --gff3 of the reference.
Steps  :
  1. keep isolates with median depth >= --min-depth.
  2. normalise each window ratio by the median ratio of its sequence (sequences >= 200 kb), which removes whole-
     sequence copy-number changes.
  3. class per window and isolate: LOW < 0.25, HALF 0.35-0.65, FULL 0.8-1.25, OTHER otherwise.
  4. candidate window: >= --min-full isolates FULL and >= --min-var isolates LOW or HALF.
     Adjacent candidate windows on one sequence merge into a region; regions >= --min-len bp are kept.
  5. null: each isolate's classes on each sequence are rotated by a random offset (keeps run lengths, breaks the
     co-location across isolates); the region count is recomputed --perms times.
  6. region annotation: share of the region covered by the other assembly in the PAF, counts of LOW/HALF/FULL by
     isolate group, median het SNPs per window in HALF and FULL isolates, and genes from the GFF3 (product / Note
     keywords HMG, homeodomain, mating, sex, pheromone).
Output : idiomorph_regions_<ref>.tsv, idiomorph_genes_<ref>.tsv, idiomorph_null_<ref>.tsv (in --out).
Limits : a lineage-private deletion also gives LOW in one lineage and FULL in the other; HALF isolates and the group
counts are given to tell them apart. Noise at 10-20x depth is why isolates below --min-depth are dropped.
"""
import argparse
import bisect
import collections
import csv
import gzip
import os
import random
import re
import statistics

LOW, HALF, FULL, OTHER = 0, 1, 2, 3
KEYWORDS = re.compile(r"HMG|homeo|mating|sex|pheromone|MAT1|SXI|STE|alpha", re.I)


def classify(r):
    return LOW if r < .25 else HALF if .35 <= r <= .65 else FULL if .8 <= r <= 1.25 else OTHER


def load_windows(path, keep):
    """-> {run: {seq: [(start, norm_ratio, het)]}}; norm by sequence median for sequences >= 200 kb."""
    raw = collections.defaultdict(lambda: collections.defaultdict(list))
    with gzip.open(path, "rt") as fh:
        next(fh)
        for line in fh:
            run, seq, start, r, snps, het = line.rstrip("\n").split("\t")
            if run in keep:
                raw[run][seq].append((int(start), float(r), int(het)))
    out = {}
    for run, d in raw.items():
        out[run] = {}
        for seq, rows in d.items():
            rows.sort()
            med = statistics.median(x[1] for x in rows) if len(rows) >= 20 else 1.0
            med = med if med > 0 else 1.0
            out[run][seq] = [(s, r / med, h) for s, r, h in rows]
    return out


def candidate_flags(classes, min_full, min_var):
    """classes: list over isolates of per-window class lists (same length) -> list of bool per window."""
    n = len(classes[0])
    flags = []
    for i in range(n):
        full = var = 0
        for c in classes:
            v = c[i]
            full += v == FULL
            var += v in (LOW, HALF)
        flags.append(full >= min_full and var >= min_var)
    return flags


def merge_regions(flags, starts, window, min_windows):
    regs, i, n = [], 0, len(flags)
    while i < n:
        if flags[i]:
            j = i
            while j + 1 < n and flags[j + 1] and starts[j + 1] == starts[j] + window:
                j += 1
            if j - i + 1 >= min_windows:
                regs.append((starts[i], starts[j] + window))
            i = j + 1
        else:
            i += 1
    return regs


def paf_intervals(path, side, min_block):
    """Intervals covered by primary alignments. side 'target' (CLFT044) or 'query' (JEL423)."""
    iv = collections.defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 12 or int(f[10]) < min_block:
                continue
            if side == "target":
                iv[f[5]].append((int(f[7]), int(f[8])))
            else:
                iv[f[0]].append((int(f[2]), int(f[3])))
    out = {}
    for k, v in iv.items():
        v.sort()
        merged = []
        for s, e in v:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[k] = merged
    return out


def covered(iv, seq, s, e):
    tot = 0
    for a, b in iv.get(seq, []):
        if b <= s:
            continue
        if a >= e:
            break
        tot += min(b, e) - max(a, s)
    return tot / (e - s)


def read_genes(gff):
    genes = collections.defaultdict(list)
    op = gzip.open if gff.endswith(".gz") else open
    with op(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
            text = " ".join(at.get(k, "") for k in ("product", "Note", "Name", "note"))
            genes[f[0]].append((int(f[3]), int(f[4]), at.get("ID", ""), text))
    for v in genes.values():
        v.sort()
    return genes


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True)
    ap.add_argument("--ref", required=True, choices=["CLFT044", "JEL423"])
    ap.add_argument("--paf", required=True)
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--groups", help="TSV run<TAB>class (isolate_classes.tsv) for group counts")
    ap.add_argument("--window", type=int, default=10000)
    ap.add_argument("--min-depth", type=float, default=20.0)
    ap.add_argument("--min-full", type=int, default=5)
    ap.add_argument("--min-var", type=int, default=3)
    ap.add_argument("--min-len", type=int, default=20000)
    ap.add_argument("--perms", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    summ = {r["run"]: r for r in csv.DictReader(open(f"{a.base}/bd_isolate_summary_{a.ref}.tsv"), delimiter="\t")}
    keep = {k for k, r in summ.items() if float(r["median_depth"]) >= a.min_depth}
    grp = {}
    if a.groups:
        for r in csv.DictReader(open(a.groups), delimiter="\t"):
            grp[r["run"]] = r["class"]
    W = load_windows(f"{a.base}/bd_windows_{a.ref}.tsv.gz", keep)
    runs = sorted(W)
    print(len(runs), "isolates with median depth >=", a.min_depth)
    seqs = sorted({s for r in runs for s in W[r]})
    minw = a.min_len // a.window
    obs, store = [], {}
    for seq in seqs:
        rs = [r for r in runs if seq in W[r] and len(W[r][seq]) == len(W[runs[0]].get(seq, W[r][seq]))]
        if len(rs) < a.min_full + a.min_var:
            continue
        starts = [x[0] for x in W[rs[0]][seq]]
        cls = [[classify(x[1]) for x in W[r][seq]] for r in rs]
        store[seq] = (rs, starts, cls)
        flags = candidate_flags(cls, a.min_full, a.min_var)
        for s, e in merge_regions(flags, starts, a.window, minw):
            obs.append((seq, s, e))
    print(len(obs), "candidate regions >=", a.min_len, "bp")
    rng = random.Random(a.seed)
    null = []
    for _ in range(a.perms):
        n = 0
        for seq, (rs, starts, cls) in store.items():
            L = len(starts)
            rot = []
            for c in cls:
                k = rng.randrange(L)
                rot.append(c[k:] + c[:k])
            n += len(merge_regions(candidate_flags(rot, a.min_full, a.min_var), starts, a.window, minw))
        null.append(n)
    with open(f"{a.out}/idiomorph_null_{a.ref}.tsv", "w") as fo:
        fo.write("perm\tregions\n")
        for i, n in enumerate(null):
            fo.write(f"{i}\t{n}\n")
    print("observed", len(obs), "null mean", round(statistics.mean(null), 2), "null max", max(null))

    iv = paf_intervals(a.paf, "target" if a.ref == "CLFT044" else "query", 5000)
    genes = read_genes(a.gff3)
    gcols = ["seq", "start", "end", "length", "windows", "aligned_to_other_ref", "n_isolates", "n_full", "n_half", "n_low",
             "n_other", "full_by_group", "half_by_group", "low_by_group", "het_median_half", "het_median_full",
             "n_genes", "keyword_genes"]
    with open(f"{a.out}/idiomorph_regions_{a.ref}.tsv", "w") as fo, open(f"{a.out}/idiomorph_genes_{a.ref}.tsv", "w") as fg:
        fo.write("\t".join(gcols) + "\n")
        fg.write("region\tgene\tstart\tend\ttext\n")
        for seq, s, e in obs:
            rs, starts, cls = store[seq]
            idx = [i for i, st in enumerate(starts) if s <= st < e]
            cnt = collections.Counter()
            bygrp = {LOW: collections.Counter(), HALF: collections.Counter(), FULL: collections.Counter()}
            hh, hf = [], []
            for ri, r in enumerate(rs):
                cs = collections.Counter(cls[ri][i] for i in idx)
                dom = max(cs, key=cs.get)  # dominant class over the region
                cnt[dom] += 1
                if dom in bygrp:
                    bygrp[dom][grp.get(r, "NA")] += 1
                het = [W[r][seq][i][2] for i in idx]
                (hh if dom == HALF else hf if dom == FULL else []).append(statistics.median(het))
            gl = [g for g in genes.get(seq, []) if g[1] > s and g[0] <= e]
            kw = [g for g in gl if KEYWORDS.search(g[3])]
            fmt = lambda c: ";".join(f"{k}:{v}" for k, v in sorted(c.items()))
            fo.write("\t".join(map(str, [seq, s, e, e - s, len(idx), round(covered(iv, seq, s, e), 3), len(rs), cnt[FULL],
                                         cnt[HALF], cnt[LOW], cnt[OTHER], fmt(bygrp[FULL]), fmt(bygrp[HALF]),
                                         fmt(bygrp[LOW]), round(statistics.median(hh), 1) if hh else "NA",
                                         round(statistics.median(hf), 1) if hf else "NA", len(gl),
                                         ";".join(g[2] for g in kw)])) + "\n")
            for g in gl:
                fg.write(f"{seq}:{s}-{e}\t{g[2]}\t{g[0]}\t{g[1]}\t{g[3][:150]}\n")


if __name__ == "__main__":
    main()
