"""Characterize default-tblastn best hits (e<=1e-5) that tblastn-fast does not reproduce.
Per lost hit: query role/gene, evalue, bitscore, pident, qcovhsp, and whether fast
found the same locus via ANY query of the same gene name, or the query at all elsewhere."""
import re, sys
from collections import Counter, defaultdict
from pathlib import Path
sys.argv = ["x"]
exec(open("/bigdata/stajichlab/jstajich/projects/MATPredict/analysis/tblastn_fast/analyze.py").read().split("timing = defaultdict")[0])

roles = {}
for f in (ROOT / "db").glob("*/**/proteins.faa"):
    rec = f.parent.name
    for l in open(f):
        m = re.match(r">([^|]+)\|gene_index=(\d+)\|name=([^|]+)\|role=(.+)$", l.strip())
        if m: roles[(m[1], m[2])] = m[4]
def role(q):
    m = re.match(r"([^|]+)\|gene(\d+)\|", q); return roles.get((m[1], m[2]), "?")

asms = sorted({p.name.split(".tblastn")[0] for p in OUT.glob("*.tblastn.tsv.zst")})
rows = []
for a in asms:
    d = load(OUT / f"{a}.tblastn.tsv.zst"); f = load(OUT / f"{a}.tblastn-fast.tsv.zst")
    best = {}
    for r in d:
        if r["qseqid"] not in best or r["bitscore"] > best[r["qseqid"]]["bitscore"]: best[r["qseqid"]] = r
    byq = defaultdict(list); byg = defaultdict(list)
    for r in f: byq[r["qseqid"]].append(r); byg[r["gene"]].append(r)
    for q, b in best.items():
        if b["evalue"] > 1e-5: continue
        L = b["e"] - b["s"] + 1
        same = lambda r: r["sseqid"] == b["sseqid"] and ov(r["s"], r["e"], b["s"], b["e"]) >= 0.5 * L
        if any(same(r) for r in byq[q]): continue
        rows.append(dict(asm=a, q=q, gene=b["gene"], role=role(q), e=b["evalue"], bits=b["bitscore"],
            pid=float(b["pident"]), cov=float(b["qcovhsp"]), alen=int(b["length"]),
            locus_by_gene=any(same(r) for r in byg[b["gene"]]), fast_has_q=bool(byq[q])))
n = len(rows); print("lost:", n)
def pct(x): return f"{x}/{n} ({100*x/n:.0f}%)"
print("by role:", Counter(r["role"] for r in rows).most_common())
print("locus still found by some fast hit of same gene name:", pct(sum(r["locus_by_gene"] for r in rows)))
print("fast has no HSP for the query anywhere:", pct(sum(not r["fast_has_q"] for r in rows)))
import statistics as st
for k in ("bits", "pid", "cov", "alen"):
    v = sorted(r[k] for r in rows); print(k, "min/median/max", v[0], st.median(v), v[-1])
bins = [(0,30),(30,39),(39,50),(50,80),(80,1e9)]
print("bitscore bins:", [(lo, sum(lo<=r["bits"]<hi for r in rows)) for lo,hi in bins])
print("evalue<=1e-10:", pct(sum(r["e"]<=1e-10 for r in rows)))
print("bits>=39 (flank floor):", pct(sum(r["bits"]>=39 for r in rows)))
print("top genes:", Counter(r["gene"] for r in rows).most_common(12))
print("per-asm lost:", Counter(r["asm"][:15] for r in rows).most_common(4))
import csv
with open(ROOT/"analysis/tblastn_fast/lost_hits.tsv","w") as fh:
    w = csv.DictWriter(fh, rows[0].keys(), delimiter="\t"); w.writeheader(); w.writerows(rows)
