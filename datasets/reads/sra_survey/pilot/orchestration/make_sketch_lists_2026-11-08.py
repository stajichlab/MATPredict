#!/usr/bin/env python3
"""Run lists for the sketch pass (2026-11-08, clone correction): every run a wave typed (any stage's runs_all.tsv, with its panel columns) that has no
<run>.sig.zip yet in the wave's out/, sketch_out/ or clone folders. Writes WAVES/<wave>/runs_sketch.tsv; prints counts.

    make_sketch_lists.py WAVES_DIR wave1 wave2 ...
"""
import csv, glob, os, sys
W = sys.argv[1]
for wave in sys.argv[2:]:
    d = os.path.join(W, wave)
    have = set()
    for sd in ("out", "sketch_out", "clone", "clone_pilot", "clone_split", "botrytis_full"):
        for p in glob.glob(os.path.join(d, sd, "**", "*.sig.zip"), recursive=True):
            have.add(os.path.basename(p)[:-8])
    rows, cols = {}, []
    for ra in sorted(glob.glob(os.path.join(d, "*", "runs_all.tsv"))):
        for r in csv.DictReader(open(ra, newline=""), delimiter="\t"):
            if not r.get("panel_dir") and not r.get("panel_prefix"): continue
            for c in r:
                if c not in cols: cols.append(c)
            rows.setdefault(r["run_accession"], r)
    todo = [r for run, r in rows.items() if run not in have]
    if todo:
        with open(os.path.join(d, "runs_sketch.tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(todo)
    print(wave, "typed runs", len(rows), "sketched", len(have & set(rows)), "to sketch", len(todo))
