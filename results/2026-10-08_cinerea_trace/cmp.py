#!/usr/bin/env python3
"""Compare two detect run folders genome by genome: every locus (called and withheld) as (family, contig, start, end, reason, genes).
usage: cmp.py RUN_A RUN_B   -> counts of identical / span-changed / reason-changed / gene-set-changed / only-in-A / only-in-B"""
import glob, os, sys, yaml
from collections import Counter

def loci(path):
    y = yaml.safe_load(open(path))
    out = []
    for x in y.get("detected") or []:
        out.append((x["family"], x["contig"], x["start"], x["end"], "called", tuple(sorted(x.get("genes_found") or []))))
    for x in y.get("suppressed_loci") or []:
        out.append((x["family"], x["contig"], x["start"], x["end"], "withheld:" + str(x.get("withheld_reason")), tuple(sorted(x.get("genes_found") or []))))
    return out

a_dir, b_dir = sys.argv[1], sys.argv[2]
tot = Counter()
for ra in sorted(glob.glob(a_dir + "/runs/*/detection_report.yaml")):
    g = ra.split("/")[-2]; rb = os.path.join(b_dir, "runs", g, "detection_report.yaml")
    if not os.path.exists(rb):
        print(g, "missing in B"); continue
    A, B = loci(ra), loci(rb)
    c = Counter()
    ka = {(l[0], l[1]): [] for l in A}; kb = {(l[0], l[1]): [] for l in B}
    for l in A: ka[(l[0], l[1])].append(l)
    for l in B: kb[(l[0], l[1])].append(l)
    for k in sorted(set(ka) | set(kb)):
        la, lb = sorted(ka.get(k, [])), sorted(kb.get(k, []))
        if la == lb: c["identical"] += len(la); continue
        if not la: c["only_in_B"] += len(lb); continue
        if not lb: c["only_in_A"] += len(la); continue
        for x, y in zip(la, lb):
            if x == y: c["identical"] += 1
            else:
                if x[2:4] != y[2:4]: c["span_changed"] += 1
                if x[4] != y[4]: c["reason_changed"] += 1
                if x[5] != y[5]: c["gene_set_changed"] += 1
        if len(la) != len(lb): c["count_changed_at_same_contig_family"] += abs(len(la) - len(lb))
    print(g, dict(c)); tot.update(c)
print("TOTAL", dict(tot))
