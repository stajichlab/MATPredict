#!/usr/bin/env python3
"""Merge per-run `matpredict reads-type` results into one table with a triage column, and list what is still missing.

Each worker writes one small TSV per run and panel (`<run>.<panel>.tsv`, see fetch_type.slurm or pilot/reads.slurm), so work can be
spread over any number of jobs and merged at any time. Nothing here recomputes a call.

triage:
  typed          one idiomorph, depth well above the noise, no warning flag (trace_* flags are noted, not blocking)
  both           both idiomorphs present (homothallic, diploid/heterokaryotic or mixed sample) - inspect depths in `flags`
  review         call is none or low_depth: candidates for assembling the hits and extending (the step we want to avoid if k-mer is enough)

Usage: accumulate.py --manifest selected.tsv --results DIR [--id-col run_accession] [--out results.tsv] [--missing missing.tsv]
"""
import argparse, csv, glob, os, sys


def triage(call, flags):
    if call == "both":
        return "both"
    if call in ("none", "low_depth", ""):
        return "review"
    return "typed"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--id-col", default="run_accession")
    ap.add_argument("--keep-cols", default="species,strain_isolate,bioproject,country,center,est_coverage,relabelled")
    ap.add_argument("--out", default="results.tsv")
    ap.add_argument("--missing", default="missing.tsv")
    a = ap.parse_args()
    meta = {r[a.id_col]: r for r in csv.DictReader(open(a.manifest), delimiter="\t")}
    keep = [c for c in a.keep_cols.split(",") if c in next(iter(meta.values()))]
    rows, seen = [], set()
    for f in sorted(glob.glob(os.path.join(a.results, "*.tsv"))):
        base = os.path.basename(f)[:-4]
        panel = base.rsplit(".", 1)[1].removesuffix("_panel") if "." in base else ""
        for r in csv.DictReader(open(f), delimiter="\t"):
            m = meta.get(r["sample"], {})
            seen.add(r["sample"])
            rows.append({"sample": r["sample"], "panel": panel, **{c: m.get(c, "") for c in keep}, "call": r["call"], "triage": triage(r["call"], r["flags"]),
                         **{k: v for k, v in r.items() if k not in ("sample", "call", "flags")}, "flags": r["flags"]})
    cols = []  # union of columns: idiomorph names (and so their depth columns) differ between panels
    for r in rows:
        cols += [k for k in r if k not in cols and k != "flags"]
    cols.append("flags")
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, restval="", delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    missing = [k for k in meta if k not in seen]
    with open(a.missing, "w") as fh:
        fh.write("\n".join(missing) + ("\n" if missing else ""))
    from collections import Counter
    print(f"{len(seen)} of {len(meta)} runs have results; {len(missing)} missing", file=sys.stderr)
    for k, v in sorted(Counter((r["panel"], r["triage"]) for r in rows).items()): print(*k, v, sep="\t", file=sys.stderr)


main()
