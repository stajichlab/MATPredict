"""Pre-sign-off regression check: list every existing call a change affects.

Curator ruling, J. Stajich, 2026-09-30: every new record or classifier rebuild
must attach this summary before sign-off. See docs/regression-check.md.

  # diff two finished run directories (any layout with */detection_report.yaml)
  python scripts/regression_check.py diff --baseline BASE_DIR --candidate CAND_DIR \
      --out OUT_DIR --title "records X+Y vs PR #9 abc123"

Several pairs can be diffed at once: repeat --pair NAME BASE CAND; each gets its
own sub-directory under --out, and --out/summary.md lists them.
"""
import argparse
import sys
from pathlib import Path

from MATPredict.detect.regression import (DEFAULT_SCORE_DELTA, compare_runs, load_reports,
                                          write_outputs)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("diff", help="diff baseline and candidate detect outputs")
    d.add_argument("--baseline", type=Path)
    d.add_argument("--candidate", type=Path)
    d.add_argument("--pair", nargs=3, action="append", metavar=("NAME", "BASE", "CAND"),
                   help="named baseline/candidate pair; repeatable")
    d.add_argument("--out", type=Path, required=True)
    d.add_argument("--title", default="regression check")
    d.add_argument("--score-delta", type=float, default=DEFAULT_SCORE_DELTA,
                   help=f"classifier margin shift (bits) that counts as a change "
                        f"(default {DEFAULT_SCORE_DELTA})")
    a = ap.parse_args(argv)

    pairs = [(n, Path(b), Path(c)) for n, b, c in (a.pair or [])]
    if a.baseline and a.candidate:
        pairs.append(("all", a.baseline, a.candidate))
    if not pairs:
        ap.error("give --baseline/--candidate or at least one --pair")
    index = [f"# {a.title}", "", "| panel | genomes | loci changed | summary |", "|---|---|---|---|"]
    for name, base, cand in pairs:
        for p in (base, cand):
            if not p.is_dir():
                print(f"missing directory: {p}", file=sys.stderr)
                return 2
        rows = compare_runs(load_reports(base), load_reports(cand), score_delta=a.score_delta)
        out = a.out / name if len(pairs) > 1 else a.out
        _, md = write_outputs(rows, out, title=f"{a.title} -- {name}", score_delta=a.score_delta)
        n_changed = sum(1 for r in rows if r["change_types"])
        index.append(f"| {name} | {len({r['genome'] for r in rows})} | {n_changed} | "
                     f"{md.relative_to(a.out)} |")
        print(f"{name}: {n_changed} changed loci -> {md}")
    if len(pairs) > 1:
        (a.out / "summary.md").write_text("\n".join(index) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
