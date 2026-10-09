#!/usr/bin/env python3
"""Withheld loci under the old and new database (Basidiomycota panel + T48-F, new code, floor 33): the question that started this.
Do the withheld spans that moved with the C. cinerea record change in cluster span but not in supported span?
usage: analyze_withheld.py (in this folder; needs newdb/ and olddb/)"""
import glob, yaml
from collections import Counter

def load(d):
    out = {}
    for f in glob.glob(f"{d}/runs/*/detection_report.yaml"):
        g = f.split("/")[-2]
        for s in yaml.safe_load(open(f)).get("suppressed_loci") or []:
            out.setdefault((g, s["family"], s["contig"]), []).append(s)
    return out

N, O = load("newdb"), load("olddb")
print("withheld loci: new db", sum(map(len, N.values())), "| old db", sum(map(len, O.values())))
cls = Counter(); rows = []
for k, a in N.items():
    b = O.get(k)
    if not b or len(a) != 1 or len(b) != 1: continue
    x, y = a[0], b[0]
    span = (x["start"], x["end"]) != (y["start"], y["end"])
    sx, sy = x.get("supported_span"), y.get("supported_span")
    cx, cy = x.get("core_span"), y.get("core_span")
    sup = (sx or {}).get("start"), (sx or {}).get("end"); sup0 = (sy or {}).get("start"), (sy or {}).get("end")
    core = (cx or {}).get("start"), (cx or {}).get("end"); core0 = (cy or {}).get("start"), (cy or {}).get("end")
    label = ("span moved" if span else "span same") + (", supported moved" if sup != sup0 else ", supported same") + (", core moved" if core != core0 else ", core same")
    cls[label] += 1
    if span: rows.append((label, k, (y["start"], y["end"]), (x["start"], x["end"]), sup0, sup, core0, core))
print("1:1 matched withheld loci:", sum(cls.values()))
for k, v in cls.most_common(): print(f"   {v:4d}  {k}")
moved = [r for r in rows if r[0].startswith("span moved")]
noise_only = [r for r in moved if "supported same" in r[0]]
print(f"\nwithheld loci whose cluster span moved: {len(moved)}; of those supported_span unchanged: {len(noise_only)} ({100*len(noise_only)/max(len(moved),1):.0f}%)")
for r in [r for r in moved if "supported moved" in r[0]][:8]:
    print("   supported moved:", r[1][0][:24], r[1][1], r[1][2], "span", r[2], "->", r[3], "supported", r[4], "->", r[5])
print("loci present in only one database:", len([k for k in N if k not in O]) + len([k for k in O if k not in N]))

print("\n--- key level (genome, family, contig): compare the SET of spans, so contigs with several loci are included ---")
keys = set(N) | set(O)
both = [k for k in keys if k in N and k in O]
def spans(L, f): return sorted(f(x) for x in L)
cl = lambda x: (x["start"], x["end"])
su = lambda x: ((x.get("supported_span") or {}).get("start"), (x.get("supported_span") or {}).get("end"))
co = lambda x: ((x.get("core_span") or {}).get("start"), (x.get("core_span") or {}).get("end"))
same_n = [k for k in both if len(N[k]) == len(O[k])]
diff_cluster = [k for k in same_n if spans(N[k], cl) != spans(O[k], cl)]
diff_cluster_sup_same = [k for k in diff_cluster if spans(N[k], su) == spans(O[k], su)]
diff_cluster_core_same = [k for k in diff_cluster if spans(N[k], co) == spans(O[k], co)]
print(f"keys in both databases: {len(both)}; with the same number of loci: {len(same_n)}; different number of loci: {len(both) - len(same_n)}; present in only one: {len(keys) - len(both)}")
print(f"keys (same count) whose cluster spans changed: {len(diff_cluster)}; supported spans unchanged in {len(diff_cluster_sup_same)} ({100*len(diff_cluster_sup_same)/max(len(diff_cluster),1):.0f}%); core spans unchanged in {len(diff_cluster_core_same)}")
print(f"keys (same count) whose supported spans changed although cluster spans did not: {len([k for k in same_n if spans(N[k], cl) == spans(O[k], cl) and spans(N[k], su) != spans(O[k], su)])}")
one = Counter()
for k in keys:
    if k in N and k not in O: one["only in new db"] += len(N[k])
    if k in O and k not in N: one["only in old db"] += len(O[k])
for k in both:
    if len(N[k]) != len(O[k]): one["loci in keys with a different count (new minus old)"] += len(N[k]) - len(O[k])
print("loci in only one database or in keys with a different count:", dict(one))
