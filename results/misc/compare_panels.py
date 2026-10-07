"""Compare two runs of the same genomes: calls by family, SLA2 in calls, adjacency.
usage: compare_panels.py OLD_DIR NEW_DIR"""
import yaml, glob, csv, collections, os, sys, statistics as st
from yaml import CSafeLoader as L
old_d, new_d = sys.argv[1], sys.argv[2]
meta = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
H = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-24_bar_synteny/hits"
def load(d):
    return {rp.split("/")[-2]: (yaml.load(open(rp), Loader=L) or {}) for rp in glob.glob(f"{d}/runs/*/detection_report.yaml")}
def orth(g):
    best = {}; p = f"{H}/{g}.tsv"
    if not os.path.exists(p): return {}
    for l in open(p):
        q, s, pid, ln, a, b, ev, bits = l.split("\t")
        if q == "none": continue
        k = q.split("|")[0]; bits = float(bits)
        if k not in best or bits > best[k][3]: best[k] = (s, min(int(a), int(b)), max(int(a), int(b)), bits)
    return best
def adj(x, o):
    return any(s == x["contig"] and max(0, max(a, x["start"]) - min(b, x["end"])) <= 20000 for s, a, b, _ in o.values())
old, new = load(old_d), load(new_d)
fam = collections.defaultdict(collections.Counter); changes = []
for g in sorted(new):
    f = meta.get(g, {}).get("ORDER") + "/" + meta.get(g, {}).get("FAMILY", "?")
    od, nd = old.get(g, {}).get("detected") or [], new[g].get("detected") or []
    c = fam[f]; c["n"] += 1; c["old"] += bool(od); c["new"] += bool(nd); o = orth(g)
    c["lost"] += bool(od and not nd); c["gained"] += bool(nd and not od)
    if nd:
        c["adj"] += any(adj(x, o) for x in nd); c["loci"] += len(nd)
        c["sla2"] += any("SLA2" in (x.get("genes_found") or []) for x in nd)
        c["high"] += sum(x["confidence"] == "high" for x in nd)
    so = sorted((x["contig"], x["start"] // 1000, x["idiomorph"]) for x in od)
    sn = sorted((x["contig"], x["start"] // 1000, x["idiomorph"]) for x in nd)
    if od and so != sn: changes.append((g, f, so, sn))
print(f"{'order/family':42s} {'n':>4} {'old':>4} {'new':>4} {'gain':>4} {'lost':>4} {'adj':>4} {'SLA2':>4} {'loci':>4} {'high':>4}")
tot = collections.Counter()
for f, c in sorted(fam.items(), key=lambda x: -x[1]["n"]):
    print(f"{f[:42]:42s} {c['n']:4d} {c['old']:4d} {c['new']:4d} {c['gained']:4d} {c['lost']:4d} {c['adj']:4d} {c['sla2']:4d} {c['loci']:4d} {c['high']:4d}")
    tot.update(c)
print(f"{'TOTAL':42s} {tot['n']:4d} {tot['old']:4d} {tot['new']:4d} {tot['gained']:4d} {tot['lost']:4d} {tot['adj']:4d} {tot['sla2']:4d} {tot['loci']:4d} {tot['high']:4d}")
print(f"\ngenomes whose EXISTING call changed place or idiomorph: {len(changes)}")
for ch in changes[:12]: print("  ", ch)
nw = collections.Counter((x["confidence"], x["locus_class"], x["idiomorph"]) for d in new.values() for x in (d.get("detected") or []))
print("new calls by confidence/class/idiomorph:", dict(nw))
w = [int(open(f"{new_d}/runs/{g}/wall_seconds").read()) for g in new if os.path.exists(f"{new_d}/runs/{g}/wall_seconds")]
if w: print(f"wall s: median {st.median(w)} p90 {sorted(w)[int(.9*len(w))]} max {max(w)}")
