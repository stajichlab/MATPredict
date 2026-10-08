#!/usr/bin/env python3
"""Fresh detect (main e9060da) versus the campaign report (v0.6.0-era), per genome: loci, idiomorphs, confidence, genes, wall time.
usage: compare_runs.py   (run in this folder)"""
import glob
import yaml

def sig(y):
    return sorted((x["contig"], x["start"], x["end"], x["idiomorph"], x["confidence"], x.get("locus_class"), ",".join(sorted(x.get("genes_found") or []))) for x in (y.get("detected") or []))

print("genome\tloci_campaign\tloci_fresh\tidentical_loci\tconf_same\twall_campaign_s\twall_fresh_s")
for d in sorted(glob.glob("inputs/*/")):
    n = d.rstrip("/").split("/")[-1]
    try:
        a = yaml.safe_load(open(f"{d}detection_report.yaml")); b = yaml.safe_load(open(f"runs/{n}/detection_report.yaml"))
    except FileNotFoundError:
        print(f"{n}\t-\t-\tno fresh run\t-\t-\t-"); continue
    sa, sb = sig(a), sig(b)
    wa = open(f"{d}wall_seconds").read().strip()
    tf = [l.split("\t") for l in open(f"out/{n}.timing.tsv")] if glob.glob(f"out/{n}.timing.tsv") else []
    wb = next((r[2] for r in tf if r[1] == "detect_fresh_s"), "")
    coords = [(s[0], s[1], s[2], s[3]) for s in sa] == [(s[0], s[1], s[2], s[3]) for s in sb]
    print(f"{n}\t{len(sa)}\t{len(sb)}\t{'yes' if sa == sb else ('coords/idiomorph same' if coords else 'DIFFERENT')}\t{[s[4] for s in sa] == [s[4] for s in sb]}\t{wa}\t{float(wb):.0f}" if wb else f"{n}\t{len(sa)}\t{len(sb)}\t{'yes' if sa == sb else 'DIFFERENT'}\t-\t{wa}\t-")
    if sa != sb:
        print("   campaign:", sa); print("   fresh:   ", sb)
