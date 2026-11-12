#!/usr/bin/env python3
"""Per-tier annotation agreement (#290): join FungiDB verdicts with each locus's confidence tier.

Input: a folder of MATPredict run folders (one `detection_report.yaml` per sample) and the
`annotation_verdicts.tsv` that `FungiDB_MATpredict/scripts/annotation_view.py` writes for the same runs.
A verdict row's `locus` is `<family with ':' -> '_'>_<position>`, where the position is the locus's index in the
report's whole `detected` list (not per family). Output: `<prefix>_by_tier.tsv` (genes, annotated, partial, novel per
tier) and `<prefix>_rows.tsv` (every verdict row with its tier and family). The annotation is a cross-check, never
the truth: `novel` means no annotated gene overlaps, not that the call is wrong.

Run it on the HPCC, e.g.  python scripts/annotation_by_tier.py <runs> <annotation_verdicts.tsv> <out-prefix>
"""
from __future__ import annotations

import argparse
import collections
import csv
import sys
from pathlib import Path

import yaml


def tiers_by_locus(runs: Path) -> dict[tuple[str, str], tuple[str, str, str | None]]:
    """{(sample, locus id): (confidence, family, idiomorph)} from every report under `runs`."""
    out = {}
    for report in sorted(runs.glob("*/detection_report.yaml")):
        doc = yaml.safe_load(report.read_text()) or {}
        for i, locus in enumerate(doc.get("detected") or []):
            fam = locus["family"].replace(":", "_")
            out[(report.parent.name, f"{fam}_{i}")] = (locus["confidence"], locus["family"], locus.get("idiomorph"))
    return out


def join(runs: Path, verdicts: Path) -> tuple[list[dict], dict[str, collections.Counter]]:
    tier_of = tiers_by_locus(runs)
    rows, counts = [], collections.defaultdict(collections.Counter)
    for r in csv.DictReader(verdicts.open(), delimiter="\t"):
        tier, family, _ = tier_of.get((r["sample"], r["locus"]), ("unmatched", "", None))
        counts[tier][r["verdict"]] += 1
        rows.append({**r, "tier": tier, "family": family})
    return rows, counts


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", type=Path)
    ap.add_argument("verdicts", type=Path)
    ap.add_argument("prefix")
    args = ap.parse_args(argv)
    rows, counts = join(args.runs, args.verdicts)
    with open(args.prefix + "_by_tier.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["tier", "genes", "annotated", "partial", "novel"])
        for tier in ("high", "medium", "low", "unmatched"):
            if tier in counts:
                c = counts[tier]
                w.writerow([tier, sum(c.values()), c["annotated"], c["partial"], c["novel"]])
    if rows:
        with open(args.prefix + "_rows.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    print(f"{len(rows)} verdict rows; unmatched to a report locus: {counts.get('unmatched', collections.Counter()).total()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
