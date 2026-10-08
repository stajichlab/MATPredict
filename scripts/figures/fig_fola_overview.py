"""Fola reads-only typing: depth space, panel effect, mixed-sample curves, blastx accuracy by reference distance.
usage: fig_fola_overview.py RESULTS_ROOT OUT.png  (RESULTS_ROOT = the worktree's results/ directory)"""
import csv
import random
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).parent))
from locusplot import C

R, OUT = Path(sys.argv[1]), Path(sys.argv[2])
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
v2 = rd(R / "2026-10-06_fola_reads_type/v2_concordance.tsv")
v4 = rd(R / "2026-10-06_fola_reads_type/v4_concordance.tsv")
mixes = rd(R / "2026-10-06_fola_reads_mixes/mixes_all.tsv")
sweep = rd(R / "2026-10-06_fola_reads_blastx/sweep_calls.tsv")
col = {"MAT1-1": C["MAT1-1"], "MAT1-2": C["MAT1-2"], "both": C["both"], "none": C["none"], "low_depth": C["none"], "MAT1": C["MAT1-1"], "MAT2": C["MAT1-2"]}
random.seed(1)
fig = plt.figure(figsize=(13.5, 9.6))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.3)

# a. depth space (Fola panel)
ax = fig.add_subplot(gs[0, 0])
for r in v4:
    ax.scatter(float(r["MAT1-2_depth"]), float(r["MAT1-1_depth"]), color=col[r["truth"]], s=18, alpha=0.7, lw=0)
for s, (dx, dy) in {"VSP-0947": (8, 8), "VSP-0931": (8, 14), "50a": (10, 14), "VSP-0980": (-10, 14)}.items():
    r = next(x for x in v4 if x["sample"] == s)
    ax.annotate(s, (float(r["MAT1-2_depth"]), float(r["MAT1-1_depth"])), xytext=(dx, dy), textcoords="offset points", fontsize=7, arrowprops=dict(arrowstyle="-", lw=0.4))
ax.set_xlabel("MAT1-2 unique k-mer depth"); ax.set_ylabel("MAT1-1 unique k-mer depth")
ax.set_title("a. 148 Fola strains, Fola-derived panel (colour: samtools breadth call)", loc="left", fontsize=9.5)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=col[k], label=k if k != "none" else "neither") for k in ("MAT1-1", "MAT1-2", "both", "none")], fontsize=7.5, frameon=False)

# b. panel effect: breadth of the carried idiomorph
ax = fig.add_subplot(gs[0, 1])
pos, ticks, labels = 0, [], []
for truth, key in (("MAT1-2", "MAT1-2_breadth"), ("MAT1-1", "MAT1-1_breadth")):
    for pn, data in (("GenBank\nreference", v2), ("Fola\nreference", v4)):
        vals = [float(r[key]) for r in data if r["truth"] == truth]
        ax.scatter([pos + random.uniform(-0.27, 0.27) for _ in vals], vals, s=7, alpha=0.5, lw=0, color=col[truth])
        med = sorted(vals)[len(vals) // 2]
        ax.plot([pos - 0.35, pos + 0.35], [med, med], color="k", lw=1.6)
        ax.text(pos, 1.04, f"{med:.2f}", ha="center", fontsize=7.5)
        ticks.append(pos); labels.append(pn); pos += 1
    ax.text(pos - 1.5, 1.12, f"{truth} strains: {truth} breadth", ha="center", fontsize=8)
    pos += 0.8
ax.set_xticks(ticks); ax.set_xticklabels(labels, fontsize=8); ax.set_ylim(0, 1.2); ax.set_ylabel("breadth of the carried idiomorph's unique k-mers")
ax.set_title("b. Panel from Fola genomes recovers the allele (black bar = median)", loc="left", fontsize=9.5)

# c. mixes
ax = fig.add_subplot(gs[1, 0])
pairs = sorted({r["pair"] for r in mixes})
ls = ["-", "--", ":"]
for pi, p in enumerate(pairs):
    sel = sorted((r for r in mixes if r["pair"] == p), key=lambda r: int(r["pct_MAT1-2"]))
    x = [int(r["pct_MAT1-2"]) for r in sel]
    ax.plot(x, [float(r["MAT1-1_depth"]) for r in sel], color=C["MAT1-1"], ls=ls[pi], lw=1.3, label=None)
    ax.plot(x, [float(r["MAT1-2_depth"]) for r in sel], color=C["MAT1-2"], ls=ls[pi], lw=1.3, label=None)
    for r in sel:
        ax.scatter(int(r["pct_MAT1-2"]), -2.3 - pi * 1.7, marker="s", s=30, color=col[r["call"]], edgecolor="none")
ax.set_xlabel("% of reads from the MAT1-2 strain"); ax.set_ylabel("unique k-mer depth"); ax.set_ylim(-7.5, 25)
ax.set_title("c. Two-strain mixtures (3 pairs, 8M reads each); squares = call", loc="left", fontsize=9.5)
ax.axhline(0, color="#888", lw=0.5)
ax.legend(handles=[Line2D([], [], color=C["MAT1-1"], label="MAT1-1 depth"), Line2D([], [], color=C["MAT1-2"], label="MAT1-2 depth")]
          + [Line2D([], [], color="k", ls=ls[i], label=p) for i, p in enumerate(pairs)]
          + [Line2D([], [], marker="s", ls="", color=col[k], label=f"call {k}") for k in ("MAT1-1", "both", "MAT1-2")], fontsize=6.5, frameon=False, ncol=2, loc="upper center")

# d. blastx by reference distance
ax = fig.add_subplot(gs[1, 1])
acc = defaultdict(lambda: [0, 0])
for r in sweep:
    k = (r["tier"], int(r["min_identity"]), r["split"]); acc[k][0] += 1; acc[k][1] += (r["truth"] == r["call"])
names = {"T1": "same species (T1)", "T2": "same genus (T2)", "T4": "other Pezizomycotina (T4)"}
tcol = {"T1": "#009E73", "T2": "#56B4E9", "T4": "#E69F00"}
for t in ("T1", "T2", "T4"):
    for split, ls_ in (("test", "-"), ("tune", ":")):
        xs = sorted({k[1] for k in acc if k[0] == t and k[2] == split})
        ax.plot(xs, [100 * acc[(t, x, split)][1] / acc[(t, x, split)][0] for x in xs], marker="o", ms=4, ls=ls_, color=tcol[t], lw=1.4)
ax.set_xlabel("minimum read-to-protein identity (%)"); ax.set_ylabel("agreement with samtools call (%)"); ax.set_ylim(0, 105)
ax.set_title("d. DIAMOND blastx of reads: accuracy by reference distance", loc="left", fontsize=9.5)
ax.legend(handles=[Line2D([], [], color=tcol[t], marker="o", label=names[t]) for t in tcol] + [Line2D([], [], color="k", ls="-", label="test half"), Line2D([], [], color="k", ls=":", label="tune half")],
          fontsize=7.5, frameon=False, loc="center left", bbox_to_anchor=(0.0, 0.45))
fig.savefig(OUT, dpi=150, bbox_inches="tight"); print("wrote", OUT)
