"""Per-family polish cap: cap6 vs off on identical code (f7b9773), same genomes."""
import collections, csv, glob, os, statistics as st
import yaml
from yaml import CSafeLoader as L
E = os.path.dirname(os.path.abspath(__file__))
L_DIR = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-24_pilots/lists"
panel = {}
for name in ("dothideomycetes", "pezizomycotina_uncurated", "orbiliomycetes", "dipodascomycetes_etc", "saccharomycetes_other"):
    for l in open(f"{L_DIR}/{name}.tsv"): panel[l.split("\t")[0]] = name
for l in open(f"{E}/genomes.tsv"):
    panel.setdefault(l.split("\t")[0], "saccharomyces_sample")
def load(v):
    out = {}
    for rp in glob.glob(f"{E}/{v}/runs/*/detection_report.yaml"):
        g = rp.split("/")[-2]; d = yaml.load(open(rp), Loader=L) or {}
        w = f"{os.path.dirname(rp)}/wall_seconds"
        out[g] = (d.get("detected") or [], int(open(w).read()) if os.path.exists(w) else None)
    return out
off, cap = load("off"), load("cap6")
def key(x): return (x["contig"], x["idiomorph"], x["locus_class"])
def ov(a, b): return a["contig"] == b["contig"] and not (a["end"] < b["start"] or a["start"] > b["end"])
by = collections.defaultdict(collections.Counter); walls = collections.defaultdict(lambda: ([], [])); diffs = []
for g in sorted(set(off) & set(cap)):
    p = panel.get(g, "?"); a, b = off[g][0], cap[g][0]; c = by[p]
    c["genomes"] += 1; c["calls_off"] += len(a); c["calls_cap"] += len(b)
    lost = [x for x in a if not any(ov(x, y) and x["idiomorph"] == y["idiomorph"] for y in b)]
    gained = [y for y in b if not any(ov(x, y) and x["idiomorph"] == y["idiomorph"] for x in a)]
    c["lost"] += len(lost); c["gained"] += len(gained)
    c["genomes_changed"] += bool(lost or gained)
    if lost or gained: diffs.append((p, g, [(x["contig"], x["start"], x["idiomorph"], x["confidence"]) for x in lost], [(y["contig"], y["start"], y["idiomorph"]) for y in gained]))
    if off[g][1] is not None and cap[g][1] is not None:
        walls[p][0].append(off[g][1]); walls[p][1].append(cap[g][1])
print(f"{'panel':26s} {'n':>4} {'calls off':>9} {'cap6':>5} {'lost':>5} {'gained':>6} {'genomes chg':>11} {'median s off':>12} {'cap6':>6} {'total h off':>11} {'cap6':>6}")
T = collections.Counter(); TW = ([], [])
for p, c in sorted(by.items()):
    wo, wc = walls[p]
    print(f"{p:26s} {c['genomes']:4d} {c['calls_off']:9d} {c['calls_cap']:5d} {c['lost']:5d} {c['gained']:6d} {c['genomes_changed']:11d} "
          f"{st.median(wo) if wo else float('nan'):12.0f} {st.median(wc) if wc else float('nan'):6.0f} {sum(wo)/3600:11.1f} {sum(wc)/3600:6.1f}")
    T.update(c); TW[0].extend(wo); TW[1].extend(wc)
print(f"{'TOTAL':26s} {T['genomes']:4d} {T['calls_off']:9d} {T['calls_cap']:5d} {T['lost']:5d} {T['gained']:6d} {T['genomes_changed']:11d} "
      f"{st.median(TW[0]):12.0f} {st.median(TW[1]):6.0f} {sum(TW[0])/3600:11.1f} {sum(TW[1])/3600:6.1f}")
print("\ngenomes whose calls differ (lost | gained):")
for d in diffs: print("  ", d)
