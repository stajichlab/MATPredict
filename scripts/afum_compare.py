#!/usr/bin/env python3
"""Compare, per strain: assembly truth (BLAST of idiomorph-specific regions), `detect` on the assembly, and `reads-type` on reads.

usage: afum_compare.py STRAIN_MAP.tsv ASSEMBLY_TRUTH.tsv DETECT_RUNS_DIR TYPES_OUT_DIR OUT_PREFIX
Truth rule: both regions present (>= 0.80 covered) = both; one present and the other absent (<= 0.20) = that idiomorph;
both absent = none; any 'partial' = ambiguous. A strain with ambiguous truth is excluded from agreement counts.
"""
from __future__ import annotations

import csv
import glob
import sys
from collections import Counter
from pathlib import Path

import yaml


def truth_call(r):
    a, b = r["MAT1-1_state"], r["MAT1-2_state"]
    if "partial" in (a, b):
        return "ambiguous"
    if a == "present" and b == "present":
        return "both"
    if a == "present":
        return "MAT1-1"
    if b == "present":
        return "MAT1-2"
    return "none"


def detect_call(runs, strain):
    f = Path(runs) / strain / "detection_report.yaml"
    if not f.exists():
        return "no_report", "", ""
    y = yaml.safe_load(open(f))
    det = y.get("detected") or []
    ids = sorted({x["idiomorph"] for x in det if x.get("idiomorph")})
    call = "none" if not ids else ids[0] if len(ids) == 1 else "both"
    conf = ",".join(x["confidence"] for x in det)
    cls = ",".join(str(x.get("locus_class")) for x in det)
    return call, conf, cls


def main(smap, truth_tsv, runs, types_dir, prefix):
    truth = {r["strain"]: r for r in csv.DictReader(open(truth_tsv), delimiter="\t")}
    reads = {}
    for f in glob.glob(f"{types_dir}/types_*.tsv"):
        for r in csv.DictReader(open(f), delimiter="\t"):
            reads[r["sample"]] = r
    rows = []
    for r in csv.DictReader(open(smap), delimiter="\t"):
        st, pref = r["strain"], r["read_prefix"]
        t = truth.get(st)
        tcall = truth_call(t) if t else "no_genome"
        dcall, conf, cls = detect_call(runs, st) if t else ("no_genome", "", "")
        rr = reads.get(pref)
        rcall = rr["call"] if rr else "no_result"
        rows.append(dict(strain=st, read_prefix=pref, layout=r["layout"], truth=tcall, detect=dcall, detect_conf=conf,
                         detect_class=cls, reads=rcall,
                         m11_depth=rr["MAT1-1_depth"] if rr else "", m12_depth=rr["MAT1-2_depth"] if rr else "",
                         shared_depth=rr["shared_depth"] if rr else "", flags=rr["flags"] if rr else ""))
    # genomes without reads
    for st, t in truth.items():
        if st not in {r["strain"] for r in rows}:
            dcall, conf, cls = detect_call(runs, st)
            rows.append(dict(strain=st, read_prefix="", layout="", truth=truth_call(t), detect=dcall, detect_conf=conf,
                             detect_class=cls, reads="no_reads", m11_depth="", m12_depth="", shared_depth="", flags=""))
    with open(f"{prefix}_comparison.tsv", "w") as o:
        w = csv.DictWriter(o, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    for label, col, subset in (("detect vs truth", "detect", lambda x: x["truth"] not in ("ambiguous", "no_genome")),
                               ("reads vs truth", "reads", lambda x: x["truth"] not in ("ambiguous", "no_genome") and x["reads"] not in ("no_reads", "no_result"))):
        sel = [x for x in rows if subset(x)]
        c = Counter((x["truth"], x[col]) for x in sel)
        agree = sum(v for (a, b), v in c.items() if a == b)
        print(f"== {label}: n={len(sel)} agree={agree} ({100 * agree / len(sel):.1f}%)")
        for k, v in sorted(c.items()):
            print(f"  truth={k[0]:9s} {col}={k[1]:10s} n={v}")
    print("truth counts:", Counter(x["truth"] for x in rows))


if __name__ == "__main__":
    main(*sys.argv[1:6])
