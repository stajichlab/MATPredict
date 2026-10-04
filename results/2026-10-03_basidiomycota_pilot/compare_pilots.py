#!/usr/bin/env python3
"""Summarise the 2026-10-03 v0.6.0 pilots (Basidiomycota, Ascomycota).

Per stratum (ORDER for Basidiomycota, CLASS for Ascomycota): genomes, reports,
genomes with >= 1 call, routing modes, wall-time median / p90 / max.
Basidiomycota also: call-by-call comparison with the 2026-09-26 full run
(run-ad1f865) on the same genomes -- the set of called families per genome.
Writes pilot_by_stratum.tsv (and basidio_vs_20260926.tsv) next to each pilot.
"""
import csv, os, statistics, sys
from collections import Counter, defaultdict
import yaml

R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"

def q(v, p):
    v = sorted(v); return v[min(len(v) - 1, int(round((len(v) - 1) * p)))] if v else ""

def load(pilot, key):
    meta = list(csv.DictReader(open(f"{R}/{pilot}/pilot_meta.tsv"), delimiter="\t"))
    out = []
    for m in meta:
        d = f"{R}/{pilot}/runs/{m['asmid']}"
        rep = f"{d}/detection_report.yaml"
        row = dict(asmid=m["asmid"], stratum=m[key] or "NA", order=m.get("order", ""), report=False, called=set(),
                   routing="", wall=None)
        if os.path.exists(f"{d}/wall_seconds"):
            row["wall"] = int(open(f"{d}/wall_seconds").read().strip() or 0)
        if os.path.exists(rep) and os.path.getsize(rep):
            doc = yaml.safe_load(open(rep)) or {}
            row["report"] = True
            row["called"] = {x.get("family") for x in (doc.get("detected") or [])}
            row["routing"] = str((doc.get("routing") or {}).get("mode", doc.get("routing_mode", "")))
        out.append(row)
    return out

def by_stratum(rows, path):
    g = defaultdict(list)
    for r in rows: g[r["stratum"]].append(r)
    with open(path, "w") as fo:
        w = csv.writer(fo, delimiter="\t", lineterminator="\n")
        w.writerow(["stratum", "genomes", "reports", "called", "routing", "wall_median_s", "wall_p90_s", "wall_max_s"])
        for s, rs in sorted(g.items(), key=lambda x: -len(x[1])):
            wall = [r["wall"] for r in rs if r["report"] and r["wall"] is not None]
            w.writerow([s, len(rs), sum(r["report"] for r in rs), sum(bool(r["called"]) for r in rs),
                        "; ".join(f"{k}:{n}" for k, n in Counter(r["routing"] for r in rs if r["report"]).most_common()),
                        int(statistics.median(wall)) if wall else "", q(wall, .9), max(wall) if wall else ""])
    wall = [r["wall"] for r in rows if r["report"] and r["wall"] is not None]
    print(f"{path}: {len(rows)} genomes, {sum(r['report'] for r in rows)} reports, "
          f"{sum(bool(r['called']) for r in rows)} with a call; wall median {statistics.median(wall):.0f} s, "
          f"p90 {q(wall, .9)} s, max {max(wall)} s, total {sum(wall)/3600:.1f} CPU-h")

b = load("2026-10-03_basidiomycota_pilot", "order")
by_stratum(b, f"{R}/2026-10-03_basidiomycota_pilot/pilot_by_stratum.tsv")
old_g = {r["genome"]: r for r in csv.DictReader(open(f"{R}/2026-09-26_basidiomycota_full/genomes.tsv"), delimiter="\t")}
old_l = defaultdict(set)
for r in csv.DictReader(open(f"{R}/2026-09-26_basidiomycota_full/loci.tsv"), delimiter="\t"):
    old_l[r["genome"]].add(r["family_called"])
cmp = Counter(); rows = []
for r in b:
    if r["asmid"] not in old_g:
        cmp["not_in_old_run"] += 1; continue
    o, n = old_l.get(r["asmid"], set()), r["called"]
    k = ("same" if o == n else "gained_call" if not o and n else "lost_call" if o and not n
         else "families_differ")
    cmp[k] += 1
    rows.append([r["asmid"], r["order"], k, "|".join(sorted(o)), "|".join(sorted(n)), old_g[r["asmid"]]["wall_s"], r["wall"]])
with open(f"{R}/2026-10-03_basidiomycota_pilot/basidio_vs_20260926.tsv", "w") as fo:
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(["asmid", "order", "change", "families_20260926", "families_v060", "wall_s_20260926", "wall_s_v060"])
    w.writerows(rows)
print("basidio vs 2026-09-26:", dict(cmp))
old_wall = [int(x[5]) for x in rows if x[5]]; new_wall = [x[6] for x in rows if x[6] is not None]
print(f"  same genomes: wall median old {statistics.median(old_wall):.0f} s, new {statistics.median(new_wall):.0f} s")

a = load("2026-10-03_ascomycota_pilot", "class")
by_stratum(a, f"{R}/2026-10-03_ascomycota_pilot/pilot_by_stratum.tsv")
