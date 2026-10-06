#!/usr/bin/env python3
"""Compare two Mucoromycota scans call by call: gained, lost, changed.

Usage: compare.py BEFORE_DIR AFTER_DIR  (each with runs/<asm>/detection_report.yaml)
"""
import glob
import os
import sys

import yaml


def load(d):
    out = {}
    for f in glob.glob(os.path.join(d, "runs", "*", "detection_report.yaml")):
        asm = f.split(os.sep)[-2]
        r = yaml.safe_load(open(f)) or {}
        calls = {}
        for x in r.get("detected") or []:
            clf = x.get("idiomorph_classifier") or {}
            best = max((clf.get("scores") or {}).values(), default=None)
            calls[(x["family"], x["contig"], x["start"])] = (
                x["idiomorph"], x["confidence"], clf.get("classifier_input"),
                None if best is None else round(best, 1))
        gated = [(s["contig"], s.get("start"), s.get("best_score"), s.get("mat_gene_min_score"), s.get("mat_gene_min_score_source"))
                 for s in r.get("suppressed_loci") or [] if s.get("reason") == "mat_gene_gate"
                 or s.get("withheld_reason") == "mat_gene_gate"]
        out[asm] = (calls, gated)
    return out


def main(before, after):
    b, a = load(before), load(after)
    both = sorted(set(b) & set(a))
    gained = lost = changed = same = 0
    lines = []
    for asm in both:
        cb, ca = b[asm][0], a[asm][0]
        for k in sorted(set(cb) | set(ca)):
            if k in cb and k not in ca:
                lost += 1
                lines.append(f"LOST    {asm} {k} {cb[k]}")
            elif k in ca and k not in cb:
                gained += 1
                lines.append(f"GAINED  {asm} {k} {ca[k]}")
            elif cb[k][:2] != ca[k][:2]:
                changed += 1
                lines.append(f"CHANGED {asm} {k} {cb[k]} -> {ca[k]}")
            else:
                same += 1
    called_b = sum(1 for x in both if b[x][0])
    called_a = sum(1 for x in both if a[x][0])
    gated_b = sum(len(b[x][1]) for x in both)
    gated_a = sum(len(a[x][1]) for x in both)
    print(f"genomes in both: {len(both)}; genomes called {called_b} -> {called_a}")
    print(f"calls unchanged {same}, gained {gained}, lost {lost}, changed {changed}")
    print(f"loci withheld by mat_gene_gate {gated_b} -> {gated_a}")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
