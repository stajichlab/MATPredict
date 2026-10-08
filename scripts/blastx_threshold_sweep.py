#!/usr/bin/env python3
"""Sweep identity filters for blastx read typing; tune on odd-numbered strains, test on even-numbered ones.

usage: blastx_threshold_sweep.py HITS_DIR TIERS_DIR CONCORDANCE_TSV OUT_PREFIX
Truth = `truth` column of the reads-type concordance table (samtools breadth call).
Per tier and filter: best hit per read (max bitscore); keep reads with identity >= I and aligned length >= MIN_LEN aa.
Class density = reads of the best-covered gene / mean reference length of that gene; control density = mean over
APN2 and SLA2. Present = density ratio to control >= MIN_RATIO and reads >= MIN_READS.
"""
from __future__ import annotations

import csv
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

MIN_LEN, MIN_RATIO, MIN_READS = 30, 0.10, 3


def ref_lengths(faa: Path):
    lens, name, n = defaultdict(list), None, 0
    for line in faa.read_text().splitlines():
        if line.startswith(">"):
            if name:
                lens[name].append(n)
            name, n = "|".join(line[1:].split("|")[:2]), 0
        else:
            n += len(line.strip())
    lens[name].append(n)
    return {k: st.mean(v) for k, v in lens.items()}


def densities(hit_file: Path, lens, min_id: float):
    best = {}
    for line in open(hit_file):
        q, s, pid, ln, *_rest, bs = line.rstrip("\n").split("\t")
        bs = float(bs)
        if q not in best or bs > best[q][0]:
            best[q] = (bs, "|".join(s.split("|")[:2]), float(pid), int(ln))
    counts = defaultdict(int)
    for bs, key, pid, ln in best.values():
        if pid >= min_id and ln >= MIN_LEN:
            counts[key] += 1
    dens = {k: v / lens[k] for k, v in counts.items()}
    return counts, dens


def call(counts, dens):
    ctrl = [dens.get(k, 0.0) for k in ("control|APN2", "control|SLA2")]
    cd = st.mean(ctrl)
    if cd <= 0:
        return "low_depth"
    present = []
    for cls in ("MAT1-1", "MAT1-2"):
        genes = [k for k in dens if k.startswith(cls + "|")]
        if genes:
            g = max(genes, key=lambda k: dens[k])
            if dens[g] / cd >= MIN_RATIO and counts[g] >= MIN_READS:
                present.append(cls)
    return "both" if len(present) == 2 else present[0] if present else "none"


def main(hits, tiers, conc, prefix):
    truth = {r["sample"]: r["truth"] for r in csv.DictReader(open(conc), delimiter="\t")}
    samples = sorted(truth)
    rows = []
    for tier in ("T1", "T2", "T4"):
        lens = ref_lengths(Path(tiers) / f"{tier}.faa")
        for min_id in (0, 40, 50, 60, 70, 80):
            for i, s in enumerate(samples):
                f = Path(hits) / f"{s}_{tier}.tsv"
                if not f.exists():
                    continue
                counts, dens = densities(f, lens, min_id)
                rows.append((tier, min_id, "tune" if i % 2 else "test", s, truth[s], call(counts, dens)))
    with open(f"{prefix}_calls.tsv", "w") as out:
        out.write("tier\tmin_identity\tsplit\tsample\ttruth\tcall\n")
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")
    summ = defaultdict(lambda: [0, 0])
    for tier, mi, split, s, t, c in rows:
        summ[(tier, mi, split)][0] += 1
        summ[(tier, mi, split)][1] += (t == c)
    print("tier min_id split n agree")
    for k in sorted(summ):
        n, a = summ[k]
        print(*k, n, a, f"{100 * a / n:.1f}%")


if __name__ == "__main__":
    main(*sys.argv[1:5])
