"""Compare tblastn default vs tblastn-fast outputs written by run_bench.sh.

Metrics, per genome and pooled:
  - wall time (median of 2 reps) and speedup
  - HSP count at e<=1e-5 (default detect cutoff is e<=10; both shown)
  - truth recall: of the curated gene placements on this assembly, how many
    are overlapped by an HSP from a query with the same gene name
  - top-hit retention: per query whose default best HSP has e<=1e-5, is that HSP
    (genome-wide) at the same place in the fast run (>=50% overlap)
"""
import json, subprocess, sys, statistics as st
from collections import defaultdict
from pathlib import Path
import yaml

ROOT = Path("/bigdata/stajichlab/jstajich/projects/MATPredict")
OUT = ROOT / "analysis/tblastn_fast/out"
cols = "qseqid sseqid pident length qlen qcovhsp evalue bitscore sstart send sframe".split()

def load(path):
    txt = subprocess.run(["zstd", "-dc", str(path)], capture_output=True, text=True).stdout
    rows = []
    for l in txt.splitlines():
        f = l.split("\t"); r = dict(zip(cols, f))
        for k in ("evalue", "bitscore"): r[k] = float(r[k])
        s, e = sorted((int(r["sstart"]), int(r["send"]))); r["s"], r["e"] = s, e
        r["gene"] = r["qseqid"].split("|")[-1]
        rows.append(r)
    return rows

def ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0) + 1)

truth_asm = json.load(open(ROOT / "results/truth_assemblies.json"))
cy = yaml.safe_load(open(ROOT / "testset/coordinate_truth/ascomycota.yaml"))
placements = defaultdict(list)  # asmid -> [(gene, acc, s, e)]
for r in cy["records"]:
    t = truth_asm.get(r["record_id"])
    if not t: continue
    for g in r["genes"]:
        for p in g["placements"]:
            placements[t["bfd_asmid"]].append((g["record_gene_name"], p["sequence_accession"], p["start"], p["end"]))
for a in placements: placements[a] = sorted(set(placements[a]))

timing = defaultdict(list); hsps = {}
for l in (OUT / "timing.tsv").read_text().splitlines()[1:]:
    a, m, s, n = l.split("\t"); timing[(a, m)].append(float(s)); hsps[(a, m)] = int(n)

asms = sorted({a for a, _ in timing if (OUT / f"{a}.tblastn-fast.tsv.zst").exists()})
tot = defaultdict(float)
print("asm\tt_def\tt_fast\tspeedup\thsp_def\thsp_fast\ttruth_def\ttruth_fast\ttop_kept")
for a in asms:
    d = load(OUT / f"{a}.tblastn.tsv.zst"); f = load(OUT / f"{a}.tblastn-fast.tsv.zst")
    def recall(rows, cut):
        hit = 0
        for g, acc, s, e in placements[a]:
            if any(r["sseqid"] == acc and r["gene"] == g and r["evalue"] <= cut and ov(r["s"], r["e"], s, e) for r in rows): hit += 1
        return hit
    n = len(placements[a])
    best = {}
    for r in d:
        if r["qseqid"] not in best or r["bitscore"] > best[r["qseqid"]]["bitscore"]: best[r["qseqid"]] = r
    kept = 0
    best = {q: b for q, b in best.items() if b["evalue"] <= 1e-5}
    for q, b in best.items():
        if any(r["qseqid"] == q and r["sseqid"] == b["sseqid"] and ov(r["s"], r["e"], b["s"], b["e"]) >= 0.5 * (b["e"] - b["s"] + 1) for r in f): kept += 1
    td, tf = st.median(timing[(a, "tblastn")]), st.median(timing[(a, "tblastn-fast")])
    for k, v in (("td", td), ("tf", tf), ("n", n), ("rd", recall(d, 10)), ("rf", recall(f, 10)), ("hd", len(d)), ("hf", len(f)), ("q", len(best)), ("kept", kept)): tot[k] += v
    print(f"{a}\t{td:.1f}\t{tf:.1f}\t{td/tf:.1f}x\t{len(d)}\t{len(f)}\t{recall(d,10)}/{n}\t{recall(f,10)}/{n}\t{kept}/{len(best)}")
print(f"TOTAL\t{tot['td']:.0f}\t{tot['tf']:.0f}\t{tot['td']/tot['tf']:.1f}x\t{int(tot['hd'])}\t{int(tot['hf'])}\t{int(tot['rd'])}/{int(tot['n'])}\t{int(tot['rf'])}/{int(tot['n'])}\t{int(tot['kept'])}/{int(tot['q'])}")
