#!/usr/bin/env python3
"""Integrity of the v0.7.0 definitive run: every listed genome has a report, or is in the first run's no-report set."""
import glob, os, csv, sys, collections
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
def listed(files):
    out = []
    for f in files:
        out += [l.split("\t")[0].strip() for l in open(f) if l.strip()]
    return out
camps = {
 "asco": (sorted(glob.glob(f"{R}/2026-10-26_v070_ascomycota/lists/wave_*.tsv")),
          f"{R}/2026-10-26_v070_ascomycota/wave_*/small_*/runs/*/detection_report.yaml",
          f"{R}/2026-10-23_ascomycota_cleanrun/wave_*/small_*/runs/*/detection_report.yaml"),
 "basidio": (sorted(glob.glob(f"{R}/2026-10-26_v070_basidiomycota/lists/*.tsv")),
          f"{R}/2026-10-26_v070_basidiomycota/runs/*/runs/*/detection_report.yaml",
          f"{R}/2026-10-23_basidiomycota_cleanrun/runs/*/runs/*/detection_report.yaml"),
}
for name, (lists, pat, oldpat) in camps.items():
    ids = listed(lists)
    new = {p.split("/")[-2] for p in glob.glob(pat)}
    old = {p.split("/")[-2] for p in glob.glob(oldpat)}
    empty = [p for p in glob.glob(pat) if os.path.getsize(p) == 0]
    miss = [g for g in ids if g not in new]
    miss_old = [g for g in ids if g not in old]
    print(f"[{name}] listed {len(ids)} unique {len(set(ids))}; reports {len(new)}; empty files {len(empty)}; "
          f"missing now {len(miss)}; missing in first run {len(miss_old)}; same missing set: {set(miss)==set(miss_old)}; "
          f"reports not in list {len(new-set(ids))}")
    print("  newly missing:", sorted(set(miss)-set(miss_old))[:10], " newly present:", sorted(set(miss_old)-set(miss))[:10])
    open(f"/tmp/miss_{name}.txt","w").write("\n".join(miss)+"\n")
