"""A. fumigatus: concordance, depth space, background breadth before/after the flank-SNP filter, depth balance of two-idiomorph strains.
usage: fig_afum_overview.py RESULTS_DIR OUT.png   (RESULTS_DIR = results/2026-10-07_afum_reads)"""
import csv
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from locusplot import C

R, OUT = Path(sys.argv[1]), Path(sys.argv[2])
rows = list(csv.DictReader(open(R / "comparison.tsv"), delimiter="\t"))
v1 = {r["sample"]: r for r in csv.DictReader(open(R / "reads_type_v1_unfiltered.tsv"), delimiter="\t")}
v2 = {r["sample"]: r for r in csv.DictReader(open(R / "reads_type_v2.tsv"), delimiter="\t")}
rows = [r for r in rows if r["reads"] not in ("no_reads", "no_result")]
col = {"MAT1-1": C["MAT1-1"], "MAT1-2": C["MAT1-2"], "both": C["both"], "none": C["none"], "low_depth": C["none"], "ambiguous": "#666666"}
random.seed(3)
fig = plt.figure(figsize=(13.5, 9.6))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.62)

# (a) call composition by assembly truth, for detect and reads
ax = fig.add_subplot(gs[0, 0])
uniq = {}
for r in rows:
    if r["truth"] in ("MAT1-1", "MAT1-2", "both") and r["strain"] not in uniq:
        uniq[r["strain"]] = r
order = ["MAT1-1", "MAT1-2", "both"]
cats = ["MAT1-1", "MAT1-2", "both", "none"]
xt, xl = [], []
for gi, t in enumerate(order):
    sel = [r for r in uniq.values() if r["truth"] == t]
    for mi, (m, key) in enumerate((("detect", "detect"), ("reads", "reads"))):
        c = Counter(r[key] if r[key] in cats else "none" for r in sel)
        bottom = 0
        for k in cats:
            v = c.get(k, 0) / len(sel)
            if v:
                ax.bar(gi * 3 + mi, v, bottom=bottom, color=col[k], edgecolor="white", width=0.9)
                if v > 0.08:
                    ax.text(gi * 3 + mi, bottom + v / 2, f"{c[k]}", ha="center", va="center", fontsize=7.5, color="white")
                bottom += v
        xt.append(gi * 3 + mi); xl.append(m)
    ax.text(gi * 3 + 0.5, 1.03, f"truth {t}\n(n={len(sel)})", ha="center", fontsize=8)
ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=8); ax.set_ylim(0, 1.18); ax.set_ylabel("fraction of strains")
ax.set_title("a. Call by method, grouped by assembly BLAST truth", loc="left", fontsize=9.5)
for k in cats:
    ax.bar(0, 0, color=col[k], label="call: " + k)
ax.legend(fontsize=7, loc="upper center", frameon=False, ncol=4, bbox_to_anchor=(0.5, -0.1))

# (b) depth space
ax = fig.add_subplot(gs[0, 1])
for r in rows:
    d1, d2 = float(r["m11_depth"] or 0), float(r["m12_depth"] or 0)
    t = r["truth"]
    if t == "no_genome":
        ax.scatter(d2, d1, facecolors="none", edgecolors=col.get(r["reads"], C["none"]), s=22, marker="^", lw=0.9)
    elif t == "ambiguous":
        ax.scatter(d2, d1, color="#666666", s=22, marker="D", alpha=0.8)
    else:
        ax.scatter(d2, d1, color=col.get(t, C["none"]), s=16, alpha=0.65, lw=0)
lab = {"IFM_59359": (8, 6), "IFM_61407": (10, -14), "DMC_AF100-1_3": (-75, -22), "AF100-12_2": (-30, 16), "F18149-Manchester": (14, -10), "AF100-12_5": (10, 8)}
for r in rows:
    if r["strain"] in lab:
        dx, dy = lab[r["strain"]]
        ax.annotate(r["strain"], (float(r["m12_depth"]), float(r["m11_depth"])), xytext=(dx, dy), textcoords="offset points", fontsize=6.8,
                    arrowprops=dict(arrowstyle="-", lw=0.4))
ax.set_xlabel("MAT1-2 unique k-mer depth"); ax.set_ylabel("MAT1-1 unique k-mer depth")
ax.set_title("b. Read depth of each idiomorph (one point per read set)", loc="left", fontsize=9.5)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=col[k], label=f"assembly {k}") for k in ("MAT1-1", "MAT1-2", "both")]
          + [Line2D([], [], marker="D", ls="", color="#666666", label="assembly partial"),
             Line2D([], [], marker="^", ls="", mfc="none", mec="k", label="no assembly (reads call colour)")], fontsize=7, frameon=False, loc="upper right")

# (c) background breadth before/after the filter
ax = fig.add_subplot(gs[1, 0])
groups = [("MAT1-2 strains\nMAT1-1 breadth", "MAT1-2", "MAT1-1_breadth"), ("MAT1-1 strains\nMAT1-2 breadth", "MAT1-1", "MAT1-2_breadth")]
pos = 0; ticks = []; labels = []
for gname, truth, key in groups:
    for vi, (vn, vd) in enumerate((("all unique\nk-mers", v1), ("flank-SNP\nfilter", v2))):
        vals = [float(vd[r["read_prefix"]][key]) for r in rows if r["truth"] == truth and r["read_prefix"] in vd]
        xs = [pos + random.uniform(-0.28, 0.28) for _ in vals]
        ax.scatter(xs, vals, s=6, color=C["MAT1-2"] if truth == "MAT1-2" else C["MAT1-1"], alpha=0.45, lw=0)
        med = sorted(vals)[len(vals) // 2]
        ax.plot([pos - 0.35, pos + 0.35], [med, med], color="k", lw=1.6)
        ticks.append(pos); labels.append(vn); pos += 1
    ax.text(pos - 1.5, 1.07, gname, ha="center", fontsize=7.5)
    pos += 0.8
hv = [min(float(v2[r["read_prefix"]]["MAT1-1_breadth"]), float(v2[r["read_prefix"]]["MAT1-2_breadth"])) for r in rows if r["truth"] == "both" and r["read_prefix"] in v2]
ax.scatter([pos + random.uniform(-0.2, 0.2) for _ in hv], hv, s=26, color=C["both"], edgecolor="k", lw=0.5)
ticks.append(pos); labels.append("lower of the\ntwo breadths"); ax.text(pos, 1.07, "assembly: both", ha="center", fontsize=7.5)
ax.set_xticks(ticks); ax.set_xticklabels(labels, fontsize=7); ax.set_ylim(0, 1.2); ax.set_ylabel("breadth of the absent idiomorph's unique k-mers")
ax.set_title("c. Background signal of the absent idiomorph, before and after the filter", loc="left", fontsize=9.5)
ax.text(0.01, 0.02, "black bar = median", transform=ax.transAxes, fontsize=7, va="bottom")

# (d) depth balance in two-idiomorph strains
ax = fig.add_subplot(gs[1, 1])
sel = [r for r in rows if r["truth"] == "both" or r["strain"] in ("DMC_AF100-1_3", "AF100-12_2")]
seen = set(); sel = [r for r in sel if not (r["strain"] in seen or seen.add(r["strain"]))]
sel.sort(key=lambda r: float(r["m11_depth"]) / max(float(r["m12_depth"]), 0.01))
y = list(range(len(sel)))
for i, r in enumerate(sel):
    d1, d2 = float(r["m11_depth"]), float(r["m12_depth"])
    ax.barh(i + 0.2, d1, height=0.38, color=C["MAT1-1"]); ax.barh(i - 0.2, d2, height=0.38, color=C["MAT1-2"])
    ax.text(max(d1, d2) + 0.6, i, f"{d1 / d2:.2f}" if d2 else "", va="center", fontsize=7)
ax.set_yticks(y); ax.set_yticklabels([f"{r['strain']}  [{r['truth'] if r['truth'] != 'no_genome' else 'no assembly'}]" for r in sel], fontsize=7)
ax.set_xlabel("unique k-mer depth (right: MAT1-1 / MAT1-2 ratio)")
ax.set_title("d. Strains with both idiomorphs (assembly or reads)", loc="left", fontsize=9.5)
ax.legend(handles=[Line2D([], [], marker="s", ls="", color=C["MAT1-1"], label="MAT1-1"), Line2D([], [], marker="s", ls="", color=C["MAT1-2"], label="MAT1-2")], fontsize=7, frameon=False, loc="lower right")
fig.savefig(OUT, dpi=150, bbox_inches="tight"); print("wrote", OUT)
