"""The A1163 MAT1-2-1 remnant across 304 A. fumigatus assemblies: minimap2 alignment presence along the A1163 locus, and remnant status by class.
usage: fig_afum_remnant.py RESULTS_DIR OUT.png   (RESULTS_DIR = results/2026-10-07_afum_reads)"""
import csv
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

sys.path.insert(0, str(Path(__file__).parent))
from locusplot import C, GENE_COLOR

R, OUT = Path(sys.argv[1]), Path(sys.argv[2])
summ = {r["strain"]: r for r in csv.DictReader(open(R / "remnant/remnant_summary.tsv"), delimiter="\t")}
strips = {l.split("\t")[0]: [float(v) for v in l.rstrip("\n").split("\t")[1].split(",")] for l in open(R / "remnant/locus_strip.tsv")}
tr = {r["strain"]: r for r in csv.DictReader(open(R / "assembly_truth.tsv"), delimiter="\t")}


def truth(s):
    a, b = tr[s]["MAT1-1_state"], tr[s]["MAT1-2_state"]
    if "partial" in (a, b):
        return "partial"
    return "both" if a == b == "present" else "MAT1-1" if a == "present" else "MAT1-2" if b == "present" else "none"


groups = ["MAT1-1", "both", "partial", "MAT1-2", "none"]
order = []
for g in groups:
    ss = [s for s in summ if truth(s) == g and s in strips]
    ss.sort(key=lambda s: -sum(strips[s]))
    order += ss
mat = np.array([strips[s] for s in order])
x0, x1 = 1653631, 1666890
fig = plt.figure(figsize=(13.5, 9.5))
gs = fig.add_gridspec(2, 2, height_ratios=[0.17, 1], width_ratios=[1.45, 1], hspace=0.04, wspace=0.5)
# gene track (A1163 scaffold_3 coordinates, ascending)
ax0 = fig.add_subplot(gs[0, 0])
genes = [("SLA2", 1653631, 1657740, "-"), ("MAT1-1-1", 1660381, 1661535, "+"), ("remnant", 1661652, 1662686, "-"), ("APN2", 1663604, 1665343, "+"), ("COX13", 1666150, 1666672, "-")]
for name, s, e, st in genes:
    xs, xe = (s - x0) / 100, (e - x0) / 100
    head = min(8, (xe - xs) * 0.4)
    pts = [(xs, -0.3), (xe - head, -0.3), (xe, 0), (xe - head, 0.3), (xs, 0.3)] if st == "+" else [(xe, -0.3), (xs + head, -0.3), (xs, 0), (xs + head, 0.3), (xe, 0.3)]
    ax0.add_patch(Polygon(pts, closed=True, facecolor=GENE_COLOR.get(name, "#ddd"), edgecolor="k", lw=0.6, hatch="////" if name == "remnant" else None))
    ax0.text((xs + xe) / 2, 0.42, name, ha="center", fontsize=7.5)
ax0.plot([0, (x1 - x0) / 100], [0, 0], color="k", lw=0.8, zorder=0)
ax0.set_xlim(0, (x1 - x0) / 100); ax0.set_ylim(-0.5, 0.8); ax0.axis("off")
ax0.set_title("a. 304 assemblies aligned to the A1163 locus (minimap2): dark = aligned, white = not", loc="left", fontsize=9.5)
# heatmap
ax1 = fig.add_subplot(gs[1, 0])
ax1.imshow(mat, aspect="auto", cmap="Greys", vmin=0, vmax=1, interpolation="nearest", extent=(0, mat.shape[1], len(order), 0))
pos = 0
for g in groups:
    n = sum(1 for s in order if truth(s) == g)
    if not n:
        continue
    ax1.axhline(pos, color=C["MAT1-1"] if g == "MAT1-1" else C["MAT1-2"] if g == "MAT1-2" else C["both"] if g == "both" else "#999", lw=0)
    dy = {"both": -14, "partial": 10}.get(g, 0)
    ax1.text(-3, pos + n / 2 + dy, f"{g} (n={n})", ha="right", va="center", fontsize=8)
    if pos:
        ax1.axhline(pos, color="#D55E00", lw=0.8)
    pos += n
ax1.set_xlim(0, (x1 - x0) / 100); ax1.set_yticks([]); ax1.set_xlabel("position along the A1163 locus (100-bp bins; SLA2 left, COX13 right)")
for name, s, e, st in genes:
    if name in ("MAT1-1-1", "remnant"):
        ax1.axvspan((s - x0) / 100, (e - x0) / 100, color="#56B4E9" if name == "remnant" else "#0072B2", alpha=0.12, lw=0)


# b. remnant status by class
ax2 = fig.add_subplot(gs[:, 1])
cats = [("full remnant, 117 bp from MAT1-1-1", "#0072B2"), ("full remnant, other spacing (insertion)", "#56B4E9"), ("3' part only (36% of remnant)", "#E69F00"), ("partial (other)", "#999999"), ("absent", "#DDDDDD")]
tab = {g: Counter() for g in groups}
for s in order:
    r = summ[s]; cov = float(r["remnant_cov"]); gap = r["remnant_mat111_gap"]
    if cov >= 0.9 and gap == "117":
        c = 0
    elif cov >= 0.9:
        c = 1
    elif 0.3 <= cov <= 0.4:
        c = 2
    elif cov > 0:
        c = 3
    else:
        c = 4
    tab[truth(s)][c] += 1
gl = [g for g in groups if sum(tab[g].values())]
for gi, g in enumerate(gl):
    n = sum(tab[g].values()); bottom = 0
    for ci, (cl, colr) in enumerate(cats):
        v = tab[g][ci]
        if v:
            ax2.barh(gi, v / n, left=bottom, color=colr, edgecolor="white", height=0.7)
            if v / n > 0.07:
                ax2.text(bottom + v / n / 2, gi, str(v), ha="center", va="center", fontsize=8, color="white" if ci in (0, 2) else "k")
            bottom += v / n
ax2.set_yticks(range(len(gl))); ax2.set_yticklabels([f"{g} (n={sum(tab[g].values())})" for g in gl], fontsize=8); ax2.invert_yaxis()
ax2.set_xlim(0, 1); ax2.set_xlabel("fraction of assemblies")
ax2.set_title("b. A1163 MAT1-2-1 remnant (1,035 bp) by assembly class", loc="left", fontsize=9.5)
ax2.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for _, c in cats], labels=[l for l, _ in cats], fontsize=7.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12))
fig.savefig(OUT, dpi=140, bbox_inches="tight"); print("wrote", OUT)
