#!/usr/bin/env python3
"""Append the runs of a typing stage (for example a top-up) to a committed wave's results_all_final.tsv and decide again the species that gained runs (2026-11-10).

Existing rows are never changed: a run already in the final table keeps its row, with its rescue and screen columns. (merge_screen.merge is for a screened stage that
covers the whole table: it rewrites every row's screen columns from the stage and blanks those of runs outside it.) New rows get an empty `rescue`. The species whose
runs were added are decided again from all their rows with wave_driver.gate_row, keeping their curator rulings and thin-evidence tokens (merge_screen.regate_species);
gates_final.tsv is rewritten and every changed row goes to --changes (default <wave>/append_changes.tsv). Standard library only.

    append_stage.py --wave-dir datasets/reads/sra_survey/pilot/wave4 --stage RESULTS/waves/wave4/thin_topup/summary/results_all.tsv [--changes FILE]
"""
import argparse, csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import merge_screen as ms  # noqa: E402

GCOLS = ["slug", "species", "decision", "reasons"]


def append_rows(final, stage):
    """(rows, added): `final` unchanged, then the stage runs it lacked (each with rescue='')."""
    seen = {r["run"] for r in final}
    new = []
    for s in stage:
        if s["run"] not in seen:
            new.append(dict(s, rescue="")); seen.add(s["run"])
    return final + new, new


def wr(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--wave-dir", required=True, help="a committed wave folder (results_all_final.tsv, gates_final.tsv)")
    ap.add_argument("--stage", required=True, help="the stage's summary/results_all.tsv")
    ap.add_argument("--changes", help="where to list the changed gate rows (default <wave-dir>/append_changes.tsv)")
    a = ap.parse_args(argv)
    fp, gp = os.path.join(a.wave_dir, "results_all_final.tsv"), os.path.join(a.wave_dir, "gates_final.tsv")
    final = ms.rd(fp)
    rows, new = append_rows(final, ms.rd(a.stage))
    wr(fp, rows, list(final[0].keys()))
    gained = {r["group"] for r in new}
    gates = ms.rd(gp)
    sub = [g for g in gates if g["slug"] in gained]
    before = {g["slug"]: (g["decision"], g["reasons"]) for g in sub}
    _, kind = ms.regate_species([r for r in rows if r["group"] in gained], sub)      # rewrites `sub` rows in place (they are the gates' dicts)
    changes = [dict(slug=g["slug"], species=g["species"], previous=before[g["slug"]][0], decision=g["decision"], reasons=g["reasons"])
               for g in sub if (g["decision"], g["reasons"]) != before[g["slug"]]]
    wr(gp, gates, GCOLS)
    wr(a.changes or os.path.join(a.wave_dir, "append_changes.tsv"), changes, ["slug", "species", "previous", "decision", "reasons"])
    print(f"{a.wave_dir}: {len(new)} runs added for {len(gained)} species; {len(changes)} gate rows changed")
    for c in changes:
        print(f"  {c['species']}: {c['previous']} -> {c['decision']}")


if __name__ == "__main__":
    main()
