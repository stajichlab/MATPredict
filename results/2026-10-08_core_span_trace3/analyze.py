#!/usr/bin/env python3
"""core_span over the 3-genome trace set: distribution of beyond_core_bp, and stability between the old and new database."""
import glob, statistics as st, yaml
def load(side):
    d = {}
    for f in glob.glob(f"{side}/runs/*/detection_report.yaml"):
        g = f.split("/")[-2]
        for x in yaml.safe_load(open(f)).get("detected") or []:
            d[(g, x["family"], x["contig"], x["core_span"]["start"] if x.get("core_span") else None, x["start"])] = x
    return d
def by_locus(side):
    d = {}
    for f in glob.glob(f"{side}/runs/*/detection_report.yaml"):
        g = f.split("/")[-2]
        for x in yaml.safe_load(open(f)).get("detected") or []:
            d.setdefault((g, x["family"], x["contig"]), []).append(x)
    return d
N, O = by_locus("newdb"), by_locus("olddb")
print("called loci: new db", sum(map(len, N.values())), "| old db", sum(map(len, O.values())))
for lab, D in (("new db", N), ("old db", O)):
    xs = [x for v in D.values() for x in v]
    b = [x["core_span"]["beyond_core_bp"] for x in xs if x.get("core_span")]
    print(f"{lab}: loci {len(xs)}, with core_span {len(b)}, beyond>0: {sum(1 for v in b if v>0)}, >=1 kb: {sum(1 for v in b if v>=1000)}, >=10 kb: {sum(1 for v in b if v>=10000)}, max {max(b) if b else 0}, median {st.median(b) if b else 0}")
    for x in xs:
        cs = x.get("core_span")
        if cs and cs["beyond_core_bp"] >= 1000: print("   ", x["family"].split(":")[1], x["contig"], x["start"], x["end"], "core", cs["start"], cs["end"], "beyond", cs["beyond_core_bp"])
print("\nloci present in both databases at the same contig and family (1:1 matches only):")
same_span = same_core = span_moved_core_same = span_moved_core_moved = 0
rows = []
for k, a in N.items():
    b = O.get(k)
    if b and len(a) == 1 and len(b) == 1:
        x, y = a[0], b[0]
        sp = (x["start"], x["end"]) == (y["start"], y["end"])
        cs = (x["core_span"] or {}).get("start"), (x["core_span"] or {}).get("end")
        co = (y["core_span"] or {}).get("start"), (y["core_span"] or {}).get("end")
        core_same = cs == co
        if sp and core_same: same_span += 1
        elif sp: same_core += 1; rows.append((k, "span same, core moved", x, y))
        elif core_same: span_moved_core_same += 1; rows.append((k, "span moved, core same", x, y))
        else: span_moved_core_moved += 1; rows.append((k, "span moved, core moved", x, y))
print(" span and core identical:", same_span, "| span moved, core unchanged:", span_moved_core_same, "| span moved and core moved:", span_moved_core_moved, "| span same, core moved:", same_core)
for k, kind, x, y in rows: print("  ", kind, k[0][:22], k[1].split(":")[1], k[2], "span", (y["start"], y["end"]), "->", (x["start"], x["end"]), "core", (y["core_span"] or {}).get("start"), (y["core_span"] or {}).get("end"), "->", (x["core_span"] or {}).get("start"), (x["core_span"] or {}).get("end"))
unmatched = [k for k in N if k not in O] + [k for k in O if k not in N]
print("loci present in only one database:", len(unmatched), [(k[0][:20], k[1].split(":")[1], k[2]) for k in unmatched][:6])
