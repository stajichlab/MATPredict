#!/usr/bin/env python3
"""Merge the rescue passes of a wave into one results table (wave 1, 2026-10-27).

Order of precedence per run: the species panel at k=31; if that is untyped (`none` or `low_depth`), the same panel at k=21 (`FALLBACK_K`; the k=21 files of the k-mer experiment
or the full2 stage); if still untyped, the genus panel (`panel_from_assemblies.py --species "Genus " --diverse`, k=31 with the k=21 fallback). A rescued row gets `rescue` = k21 or
genus_panel; its flags, depths and call are those of the rescuing pass. Gate the result with `wave_driver.py gate --stage final` (put the output in <wave>/final/summary/).

    merge_rescues.py --base results_all.tsv [--base2 results_all_full.tsv] [--k21-dir DIR --k21-runs kexp_runs.tsv] [--full2 results_all.tsv] [--genus results_all.tsv] --out results_all_final.tsv
"""
import argparse, csv, os

COLS = ("call", "depth_1", "depth_2", "shared_depth", "breadth_1", "breadth_2", "flags")


def rd(p):
    return list(csv.DictReader(open(p, newline="").read().replace("\r", "").splitlines(), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--base", required=True); ap.add_argument("--base2", help="results that take precedence over --base for the same run (the full stage)")
    ap.add_argument("--k21-dir"); ap.add_argument("--k21-runs"); ap.add_argument("--full2"); ap.add_argument("--genus"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    base = {r["run"]: r for r in rd(a.base)}
    if a.base2: base.update({r["run"]: r for r in rd(a.base2)})
    for r in base.values(): r["rescue"] = ""
    untyped = lambda r: r["state"] == "done" and r["call"] in ("none", "low_depth")
    ok_names = lambda r: tuple((r.get("idiomorphs") or "MAT1-1 MAT1-2").split()) + ("both",)
    if a.k21_dir and a.k21_runs:
        for e in rd(a.k21_runs):
            r = base.get(e["run_accession"]); p = os.path.join(a.k21_dir, f"{e['run_accession']}.{e['panel_prefix']}.k21.tsv")
            if r and untyped(r) and os.path.exists(p):
                x = rd(p)[0]
                if x["call"] in ok_names(r):
                    r.update(call=x["call"], depth_1=x[(r.get("idiomorphs") or "MAT1-1 MAT1-2").split()[0] + "_depth"], depth_2=x[(r.get("idiomorphs") or "MAT1-1 MAT1-2").split()[1] + "_depth"], shared_depth=x["shared_depth"], breadth_1=x[(r.get("idiomorphs") or "MAT1-1 MAT1-2").split()[0] + "_breadth"], breadth_2=x[(r.get("idiomorphs") or "MAT1-1 MAT1-2").split()[1] + "_breadth"], flags=x["flags"], rescue="k21")
    for path, tag in ((a.full2, "k21"), (a.genus, "genus_panel")):
        if path:
            for x in rd(path):
                r = base.get(x["run"])
                if r and untyped(r) and x["call"] in ok_names(r):
                    r.update({k: x[k] for k in COLS}); r["rescue"] = tag
    rows = list(base.values())
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0].keys()), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print(len(rows), "runs;", sum(1 for r in rows if r["rescue"]), "rescued")


if __name__ == "__main__":
    main()
