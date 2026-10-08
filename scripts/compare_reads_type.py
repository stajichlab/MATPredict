#!/usr/bin/env python3
"""Compare `matpredict reads-type` TSVs with samtools-breadth calls (breadth >= 90% of an idiomorph reference).

usage: compare_reads_type.py MAT_coverage.tsv OUT_PREFIX types_*.tsv
MAT_coverage.tsv columns: Isolate, MAT-1, MAT-2 (percent breadth).
"""
from __future__ import annotations

import csv
import sys
from collections import Counter


def truth_call(m1: float, m2: float) -> str:
    a, b = m1 >= 90, m2 >= 90
    return "both" if a and b else "MAT1-1" if a else "MAT1-2" if b else "none"


def main(cov: str, prefix: str, tsvs: list[str]) -> None:
    truth = {}
    with open(cov) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            truth[row["Isolate"]] = (truth_call(float(row["MAT-1"]), float(row["MAT-2"])), row["MAT-1"], row["MAT-2"])
    rows = []
    for t in tsvs:
        with open(t) as fh:
            rows += list(csv.DictReader(fh, delimiter="\t"))
    table = Counter()
    with open(f"{prefix}_concordance.tsv", "w") as out:
        out.write("sample\ttruth\tkmer_call\tagree\tsamtools_MAT1_breadth\tsamtools_MAT2_breadth\tMAT1-1_breadth\tMAT1-2_breadth\tMAT1-1_depth\tMAT1-2_depth\tshared_depth\tflags\n")
        for r in sorted(rows, key=lambda r: r["sample"]):
            tr, b1, b2 = truth[r["sample"]]
            table[(tr, r["call"])] += 1
            out.write("\t".join([r["sample"], tr, r["call"], str(tr == r["call"]), b1, b2, r["MAT1-1_breadth"], r["MAT1-2_breadth"],
                                 r["MAT1-1_depth"], r["MAT1-2_depth"], r["shared_depth"], r["flags"]]) + "\n")
    n = sum(table.values())
    agree = sum(v for (a, b), v in table.items() if a == b)
    print(f"n={n} agree={agree} ({100 * agree / n:.1f}%)")
    for (a, b), v in sorted(table.items()):
        print(f"truth={a:7s} kmer={b:10s} n={v}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
