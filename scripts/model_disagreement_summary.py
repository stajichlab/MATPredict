"""Summarise the census of `polished_disagree` genes against two reference-based yardsticks (#290, #349).

Input: one or more `<run>.census.tsv` from `model_disagreement_census.py`. For each run and subset (all genes, core_MAT genes, the
`trigger_noncanonical` genes, and each census class) it counts the genes and, over the genes with both identities and a reference length,
how often the ALTERNATE model has the higher identity to the reference, how often its length is closer to the reference length, the
median identity difference (alternate minus primary, percentage points) and the median length difference (aa). Identity and length come
from the same tools that made the models, so they are not independent of them; they are independent of the splice-motif check, which is why
they are used together.

Usage: model_disagreement_summary.py OUT.tsv CENSUS.tsv [CENSUS.tsv ...]
"""
from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def yardsticks(rows: list[dict]) -> dict:
    """Counts and medians over `rows` (census rows as dicts of str)."""
    ident = [(_f(r["alternate_identity"]), _f(r["primary_identity"])) for r in rows]
    ident = [(a, p) for a, p in ident if a is not None and p is not None]
    length = [(_f(r["alternate_aa"]), _f(r["primary_aa"]), _f(r["ref_aa"])) for r in rows]
    length = [(a, p, ref) for a, p, ref in length if None not in (a, p, ref) and ref > 0]
    out = {"genes": len(rows), "with_identity": len(ident), "with_ref_length": len(length)}
    out["alt_higher_identity"] = sum(1 for a, p in ident if a > p)
    out["identity_tie"] = sum(1 for a, p in ident if a == p)
    out["median_identity_diff"] = round(statistics.median(a - p for a, p in ident), 2) if ident else ""
    out["alt_closer_length"] = sum(1 for a, p, ref in length if abs(a - ref) < abs(p - ref))
    out["length_tie"] = sum(1 for a, p, ref in length if abs(a - ref) == abs(p - ref))
    out["median_aa_diff"] = round(statistics.median(a - p for a, p, ref in length), 1) if length else ""
    return out


COLS = ["run", "subset", "genes", "with_identity", "alt_higher_identity", "identity_tie", "median_identity_diff", "with_ref_length", "alt_closer_length", "length_tie", "median_aa_diff"]


def run(paths: list[Path], out: Path) -> list[dict]:
    rows_out = []
    for path in paths:
        label = path.name.split(".")[0]
        with path.open() as fh:
            rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["class"] != "not_checked"]
        subsets = {"all": rows, "core_MAT": [r for r in rows if r["role"] == "core_MAT"],
                   "trigger_noncanonical": [r for r in rows if r["trigger_noncanonical"] == "True"]}
        for c in ("primary_issue_alt_clean", "primary_clean_alt_issue", "both_clean", "both_issue"):
            subsets[c] = [r for r in rows if r["class"] == c]
            subsets[c + "|core_MAT"] = [r for r in subsets[c] if r["role"] == "core_MAT"]
        for name, rs in subsets.items():
            rows_out.append({"run": label, "subset": name, **yardsticks(rs)})
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows_out)
    return rows_out


def main(argv: list[str]) -> int:
    run([Path(p) for p in argv[1:]], Path(argv[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
