#!/usr/bin/env python3
"""Compare the pilot's `matpredict detect` calls (assemblies) with `reads-type` calls (reads), per strain and panel.
Usage: compare.py pilot_samples.tsv pilot_out > comparison.tsv   (pilot_out has detect/runs/<strain>/detection_report.yaml and reads/<strain>.<panel>_panel.tsv)"""
import csv, os, sys, yaml

samples, out = sys.argv[1:3]
w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
w.writerow(["strain", "species", "relabelled", "center", "detect_call", "detect_confidence", "detect_n_loci", "panel", "reads_call", "MAT1-1_depth", "MAT1-2_depth", "shared_depth", "flags", "agree"])
for s in csv.DictReader(open(samples), delimiter="\t"):
    rep = yaml.safe_load(open(f"{out}/detect/runs/{s['strain']}/detection_report.yaml"))
    det = rep.get("detected") or []
    ids = sorted({d.get("idiomorph") for d in det if d.get("idiomorph") not in (None, "undetermined")})
    dcall = "+".join(ids) if ids else ("undetermined" if det else "not_detected")
    conf = max((d.get("confidence") for d in det), key=lambda c: {"high": 3, "medium": 2, "low": 1}.get(c, 0), default="")
    for panel in ("posadasii", "immitis"):
        r = next(csv.DictReader(open(f"{out}/reads/{s['strain']}.{panel}_panel.tsv"), delimiter="\t"))
        agree = "yes" if r["call"] == dcall else ("both-vs-one" if r["call"] == "both" or "+" in dcall else "no")
        w.writerow([s["strain"], s["species"].split()[-1], s["relabelled"], s["center"], dcall, conf, len(det), panel, r["call"], r["MAT1-1_depth"], r["MAT1-2_depth"], r["shared_depth"], r["flags"], agree])
