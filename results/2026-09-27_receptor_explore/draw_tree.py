"""Named STE3 receptor tree PDF: tips = species + strain; curated mating receptors bold and coloured."""
import csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo
hits = {r["seq_id"]: r for r in csv.DictReader(open("ste3_hits.tsv"), delimiter="\t")}
COL = {"Agaricales_PR": "#c41e1e", "Ustilaginales_pra": "#1f4fbf", "Tremellales_STE3": "#6a3d9a",
       "Sporidiobolales": "#e6550d", "Wallemia_STE3": "#2ca25f", "OUT": "#000000"}
ORDERCOL = {"Agaricales": "#444444", "Pucciniales": "#b8860b", "Sporidiobolales": "#e6550d"}
def grp(n):
    if n.startswith("OUT|"): return "OUT"
    g = n.split("|")[2]
    if g in ("bar3", "bbr2", "pheromone_receptor"): return "Agaricales_PR"
    if g == "pra1": return "Ustilaginales_pra"
    if g in ("STE3a1", "STE3a2"): return "Sporidiobolales"
    if "wallMAT" in n: return "Wallemia_STE3"
    return "Tremellales_STE3"
t = Phylo.read("tree_ft.treefile", "newick")
t.root_with_outgroup([x for x in t.get_terminals() if x.name.startswith("OUT|")][0])
t.ladderize(reverse=True)
terms = t.get_terminals(); n = len(terms)
y = {c: i for i, c in enumerate(terms)}
def yof(c):
    if c in y: return y[c]
    v = [yof(k) for k in c.clades]; y[c] = (min(v) + max(v)) / 2; return y[c]
yof(t.root); d = t.depths(); mx = max(d.values())
fig, ax = plt.subplots(figsize=(14, max(8, n * 0.13 + 2)))
for c in t.find_clades():
    if c.clades:
        ys = [y[k] for k in c.clades]; ax.plot([d[c], d[c]], [min(ys), max(ys)], "k-", lw=0.4)
        for k in c.clades: ax.plot([d[c], d[k]], [y[k], y[k]], "k-", lw=0.4)
        if c.confidence is not None and c.confidence >= 0.9 and c != t.root: ax.plot(d[c], y[c], "ko", ms=1.6)
for c in terms:
    if c.name.startswith(("REF|", "OUT|")):
        g = grp(c.name); lab = ("S. cerevisiae STE3 (root)" if g == "OUT" else
                                f"REF {c.name.split('|')[2]}  {c.name.split('|')[1]}")
        ax.text(d[c] + mx * 0.004, y[c], lab, fontsize=5.2, va="center", color=COL[g], fontweight="bold")
    else:
        h = hits[c.name]
        ax.text(d[c] + mx * 0.004, y[c], f"{h['species']} {h['strain']}  [{h['order']}]".strip(), fontsize=5.0,
                va="center", color=ORDERCOL.get(h["order"], "#777777"))
ax.set_ylim(n, -1); ax.set_xlim(0, mx * 1.5); ax.axis("off")
ax.plot([0, 0.5], [n - 0.2, n - 0.2], "k-", lw=0.8); ax.text(0.25, n + 0.8, "0.5 subst/site", ha="center", fontsize=6)
hd = [Line2D([], [], color=v, marker="s", lw=0, ms=6, label=k) for k, v in COL.items()]
hd += [Line2D([], [], color=ORDERCOL["Pucciniales"], marker="s", lw=0, ms=6, label="query: Pucciniales"),
       Line2D([], [], color="k", marker="o", lw=0, ms=3, label="FastTree local support >= 0.9")]
ax.legend(handles=hd, loc="upper left", fontsize=7, frameon=False)
ax.set_title("Basidiomycota STE3-like receptors (Pfam PF02076, hmmalign 296 cols, FastTree -lg -gamma); curated mating receptors in bold", fontsize=9, loc="left")
fig.savefig("receptor_tree.pdf", bbox_inches="tight"); print("tips", n)
