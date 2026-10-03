"""Split BFD Basidiomycota into `short`-partition jobs of ~1.25 h real runtime.

Per-genome cost comes from the gate C pilot (results/2026-09-26_basidio_gateC):
median wall for lineage-routed genomes (orders inside a curated family scope)
and for phylum_fallback genomes (every other order). A job runs SLOTS genomes
at once (xargs -P = cpus), so its expected wall is sum(cost)/SLOTS. Genomes are
packed per order, largest file first, so the slow genomes start early.

usage: make_chunks.py IN_SCOPE_S FALLBACK_S [SLOTS] [TARGET_S]
"""
import csv, os, sys
IN_S, FB_S = float(sys.argv[1]), float(sys.argv[2])
SLOTS = int(sys.argv[3]) if len(sys.argv) > 3 else 16
TARGET = float(sys.argv[4]) if len(sys.argv) > 4 else 4500   # 1.25 h
SCOPED = {"Agaricales", "Russulales", "Boletales", "Polyporales", "Tremellales",
          "Ustilaginales", "Malasseziales"}                  # run-ad1f865 db scopes
# Agaricales is searched against five families (HD, PR, Aalpha, Balpha, Bbeta):
# pilot walls 209, 263, 71, 56 s (mean ~150 s), not the 68 s lineage mean.
ORDER_COST = {"Agaricales": 160.0}
S = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"
L = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lists")
by_order = {}
for r in csv.DictReader(open(S)):
    if r["PHYLUM"] != "Basidiomycota" or not r["NCBI_TAXONID"].strip():
        continue
    f = os.path.join(L, r["ASMID"] + ".fa.gz")
    if not os.path.exists(f):
        continue
    by_order.setdefault(r["ORDER"] or "unassigned", []).append((os.path.getsize(f), r["ASMID"], r["NCBI_TAXONID"].strip()))
chunks = []           # (name, [(asmid, taxid)], expected_s)
small = []            # orders too small for their own job are pooled
for order, gs in sorted(by_order.items(), key=lambda kv: -len(kv[1])):
    cost = ORDER_COST.get(order, IN_S if order in SCOPED else FB_S)
    per_job = max(1, int(TARGET * SLOTS / cost))
    gs.sort(reverse=True)
    if len(gs) * cost / SLOTS < TARGET / 3:
        small.extend((cost, g) for g in gs); continue
    n = -(-len(gs) // per_job)
    for i in range(n):
        part = gs[i::n]   # interleave so each chunk gets a share of the big genomes
        chunks.append((f"{order}_{i+1}of{n}", [(a, t) for _, a, t in part], len(part) * cost / SLOTS))
small.sort(key=lambda x: -x[1][0])
cur, cur_s, k = [], 0.0, 1
for cost, (_, a, t) in small:
    if cur and cur_s + cost / SLOTS > TARGET:
        chunks.append((f"pooled_{k}", cur, cur_s)); cur, cur_s, k = [], 0.0, k + 1
    cur.append((a, t)); cur_s += cost / SLOTS
if cur: chunks.append((f"pooled_{k}", cur, cur_s))
os.makedirs(OUT, exist_ok=True)
tot = 0
for name, gs, s in chunks:
    safe = "".join(c if c.isalnum() or c in "_-" else "_" for c in name)
    with open(os.path.join(OUT, safe + ".tsv"), "w") as fh:
        for a, t in gs: fh.write(f"{a}\t{t}\n")
    tot += len(gs)
    print(f"{safe}\t{len(gs)}\t{s/3600:.2f}h")
print(f"TOTAL\t{tot} genomes in {len(chunks)} jobs", file=sys.stderr)
