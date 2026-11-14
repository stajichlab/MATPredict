#!/usr/bin/env python3
"""Isolate-specific chromosome depth changes (possible aneuploidy or copy-number change) in Bd, from bd_depth_by_seq_<ref>.tsv.

Input  --depth-by-seq   bd_depth_by_seq_<ref>.tsv (run, seq, length, depth_ratio, het_per_mb; ratio = sequence median depth / genome median)
       --classes        tables/isolate_classes.tsv (run, isolate, class, ...; class from the parent-of-origin analysis)
       --summary        bd_isolate_summary_<ref>.tsv (median_depth per run)
Method : sequences of at least --min-len bp. For each class with at least --min-class isolates, the median ratio of a sequence over
         that class is the reference level (this removes reference-specific offsets, for example a lineage that maps worse on one
         chromosome). Normalised ratio = ratio / class median. An isolate-sequence is flagged when the normalised ratio is outside
         1 +/- --delta. Only isolates with median depth >= --min-depth are flagged (low depth is noisy). Expected values for one extra
         or one missing copy in a diploid are 1.5 and 0.5; values near 1.33/0.67 fit a 3n or 1n state of a 4n/2n genome.
With --fallback-class the isolates of small classes are compared with the median of that class (labelled with a star).
Output : aneuploidy_flags.tsv (run, isolate, class, seq, length, ratio, class_median, normalised, het_per_mb, median_depth) and
         aneuploidy_by_sequence.tsv (seq, length, per class median ratio and number of flagged isolates).
Limit  : depth ratio is relative to the genome median, so a change of one chromosome shifts the others slightly; a whole-genome ploidy
         change is invisible here. Repeat-rich or organellar sequences (ratio far from 1 in every isolate) are excluded by --min-len and
         by --max-class-median-dev.
"""
import argparse
import collections
import csv
import statistics


def normalise(ratios_by_class):
    """{class: {run: ratio}} -> ({class: median}, {(class, run): ratio / median})."""
    med, norm = {}, {}
    for c, d in ratios_by_class.items():
        m = statistics.median(d.values())
        med[c] = m
        for run, r in d.items():
            norm[(c, run)] = r / m if m else float("nan")
    return med, norm


def flag(norm_value, delta):
    return abs(norm_value - 1.0) >= delta


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--depth-by-seq", required=True)
    ap.add_argument("--classes", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-len", type=int, default=500000)
    ap.add_argument("--delta", type=float, default=0.25)
    ap.add_argument("--min-class", type=int, default=10)
    ap.add_argument("--min-depth", type=float, default=20.0)
    ap.add_argument("--max-class-median-dev", type=float, default=0.5)
    ap.add_argument("--fallback-class", help="class whose median is the reference level for isolates of classes with fewer than --min-class isolates (for example JEL-type for the two HET-rich isolates); such flags are labelled with the class name plus '*'")
    a = ap.parse_args()
    cls = {r["run"]: r for r in csv.DictReader(open(a.classes), delimiter="\t")}
    depth = {r["run"]: float(r["median_depth"]) for r in csv.DictReader(open(a.summary), delimiter="\t")}
    rows = [r for r in csv.DictReader(open(a.depth_by_seq), delimiter="\t") if int(r["length"]) >= a.min_len and r["run"] in cls]
    ncls = collections.Counter(cls[r["run"]]["class"] for r in {x["run"]: x for x in rows}.values())
    usable = {c for c, n in ncls.items() if n >= a.min_class}
    seqs = sorted({r["seq"] for r in rows})
    flags, bysec = [], []
    for s in seqs:
        per = collections.defaultdict(dict)
        info = {}
        for r in rows:
            if r["seq"] == s and cls[r["run"]]["class"] in usable:
                per[cls[r["run"]]["class"]][r["run"]] = float(r["depth_ratio"])
                info[r["run"]] = (r["length"], r["het_per_mb"])
        med, norm = normalise(per)
        if a.fallback_class and a.fallback_class in med:
            small = collections.defaultdict(dict)
            for r in rows:
                c = cls[r["run"]]["class"]
                if r["seq"] == s and c not in usable:
                    small[c + "*"][r["run"]] = float(r["depth_ratio"])
                    info[r["run"]] = (r["length"], r["het_per_mb"])
            for c, d in small.items():
                med[c] = med[a.fallback_class]
                for run, v in d.items():
                    norm[(c, run)] = v / med[a.fallback_class]
                    per[c] = d
        nflag = collections.Counter()
        for (c, run), v in norm.items():
            if abs(med[c] - 1.0) > a.max_class_median_dev:
                continue
            if flag(v, a.delta) and depth.get(run, 0) >= a.min_depth:
                nflag[c] += 1
                flags.append((run, cls[run]["isolate"], c, s, info[run][0], round(per[c][run], 3), round(med[c], 3), round(v, 3),
                              info[run][1], depth[run]))
        bysec.append((s, rows[[r["seq"] for r in rows].index(s)]["length"],
                      ";".join(f"{c}:{med[c]:.2f}:{nflag[c]}" for c in sorted(med))))
    with open(a.out + "/aneuploidy_flags.tsv", "w") as fo:
        fo.write("run\tisolate\tclass\tseq\tlength\tratio\tclass_median\tnormalised\thet_per_mb\tmedian_depth\n")
        for f in sorted(flags, key=lambda x: (x[0], x[3])):
            fo.write("\t".join(map(str, f)) + "\n")
    with open(a.out + "/aneuploidy_by_sequence.tsv", "w") as fo:
        fo.write("seq\tlength\tclass:median_ratio:n_flagged\n")
        for b in bysec:
            fo.write("\t".join(b) + "\n")
    print(len(flags), "isolate-sequence flags;", len({f[0] for f in flags}), "isolates; classes used:", sorted(usable))


if __name__ == "__main__":
    main()
