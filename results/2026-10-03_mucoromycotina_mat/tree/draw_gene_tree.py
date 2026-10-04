#!/usr/bin/env python3
"""Draw and summarise a sexP/sexM gene tree (IQ-TREE .treefile, labels SH-aLRT/UFBoot).

Rooting (curator ruling 2026-10-03): non-MAT HMG-box outgroup (tips OUT_*). The tree is
rooted on the outgroup tip farthest from the ingroup, then, if the outgroup tips form
one clade, on that clade.

Reported (gene_tree_summary_<aln>.tsv):
  - ingroup (sexP + sexM) monophyletic?  support at its MRCA
  - sexP tips monophyletic? support; number of sexM tips inside the sexP MRCA, and
    the reverse
  - per tip: idiomorph label (classifier call) and whether it sits in its own clade
Figure gene_tree_<aln>.{png,svg}: rectangular tree, outgroup collapsed to one
wedge, tips coloured by call (sexP red, sexM blue), curated records marked, labels
"Genus species strain (source)", UFBoot >= 95 marked on internal nodes.
"""
import csv, re, sys
from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat/tree"
aln = sys.argv[1] if len(sys.argv) > 1 else "hmg"
tips = {r["tip"]: r for r in csv.DictReader(open(f"{T}/tips.tsv"), delimiter="\t")}
tree = Phylo.read(f"{T}/iq_{aln}.treefile", "newick")

def sup(c):
    """(SH-aLRT, UFBoot) from an IQ-TREE label 'a/b', else (None, None)."""
    lab = c.name if c.name and "/" in str(c.name) else (str(c.confidence) if c.confidence is not None else None)
    if lab and "/" in lab:
        a, b = lab.split("/")[:2]
        return float(a), float(b)
    return None, None

leaves = tree.get_terminals()
out = [t for t in leaves if t.name.startswith("OUT_")]
ing = [t for t in leaves if not t.name.startswith("OUT_")]
P = [t for t in ing if tips.get(t.name, {}).get("idiomorph") == "Plus"]
M = [t for t in ing if tips.get(t.name, {}).get("idiomorph") == "Minus"]
far = max(out, key=lambda t: tree.distance(t, ing[0]))
tree.root_with_outgroup(far)
oc = tree.common_ancestor(out)
if set(oc.get_terminals()) == set(out):
    tree.root_with_outgroup(oc)
tree.ladderize()

def clade_stats(group, other):
    m = tree.common_ancestor(group)
    inside = set(m.get_terminals())
    return dict(n=len(group), mrca_tips=len(inside), others_inside=len([t for t in other if t in inside]),
                outgroup_inside=len([t for t in out if t in inside]), shalrt=sup(m)[0], ufboot=sup(m)[1],
                monophyletic=len(inside) == len(group))
rows = []
for name, g, o in (("ingroup", ing, out), ("sexP", P, M), ("sexM", M, P)):
    s = clade_stats(g, o); s["set"] = name; rows.append(s)
with open(f"{T}/gene_tree_summary_{aln}.tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=["set", "n", "mrca_tips", "others_inside", "outgroup_inside", "monophyletic", "shalrt", "ufboot"],
                       delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)
for r in rows:
    print(r)

# figure: collapse the outgroup clade(s) to a single labelled tip
if set(oc.get_terminals()) == set(out):
    n_out = len(out)
    oc.clades = []
    oc.name = f"OUTGROUP ({n_out} non-MAT HMG-box proteins)"
leaves = tree.get_terminals()
depth = tree.depths()
if max(depth.values()) == 0:
    depth = tree.depths(unit_branch_lengths=True)
ypos = {t: i for i, t in enumerate(leaves)}
def y_of(c):
    if c.is_terminal():
        return ypos[c]
    ys = [y_of(k) for k in c.clades]
    c._y = (min(ys) + max(ys)) / 2
    return c._y
y_of(tree.root)
H = max(8, 0.075 * len(leaves))
fig, ax = plt.subplots(figsize=(10, H))
def draw(c):
    x, y = depth[c], (ypos[c] if c.is_terminal() else c._y)
    for k in c.clades:
        xk, yk = depth[k], (ypos[k] if k.is_terminal() else k._y)
        ax.plot([x, x], [y, yk], color="black", lw=0.4)
        ax.plot([x, xk], [yk, yk], color="black", lw=0.4)
        if not k.is_terminal():
            a, b = sup(k)
            if b is not None and b >= 95:
                ax.plot(xk, yk, marker="o", ms=1.6, color="black")
        draw(k)
draw(tree.root)
xmax = max(depth.values())
COL = {"Plus": "#d62728", "Minus": "#1f77b4"}
for t in leaves:
    r = tips.get(t.name)
    if r is None:
        ax.text(depth[t] + xmax * 0.01, ypos[t], t.name, va="center", fontsize=6, color="grey")
        continue
    col = COL.get(r["idiomorph"], "grey")
    lab = re.sub(r"\s+", " ", r["label"])
    weight = "bold" if r["kind"] == "record" else "normal"
    ax.plot(depth[t], ypos[t], marker="s", ms=2, color=col)
    ax.text(depth[t] + xmax * 0.01, ypos[t], lab, va="center", fontsize=3.2, color=col, fontweight=weight)
ax.set_ylim(len(leaves), -1)
ax.axis("off")
s = {r["set"]: r for r in rows}
ax.set_title(f"sexP / sexM gene tree ({'HMG box' if aln == 'hmg' else 'full length'}; IQ-TREE, UFBoot 1000)\n"
             f"sexP: {s['sexP']['n']} tips, MRCA holds {s['sexP']['others_inside']} sexM, UFBoot {s['sexP']['ufboot']}; "
             f"sexM: {s['sexM']['n']} tips, MRCA holds {s['sexM']['others_inside']} sexP, UFBoot {s['sexM']['ufboot']}",
             fontsize=8)
ax.plot([], [], "s", color=COL["Plus"], label="Plus call (sexP)")
ax.plot([], [], "s", color=COL["Minus"], label="Minus call (sexM)")
ax.plot([], [], "o", color="black", ms=3, label="UFBoot >= 95")
ax.legend(loc="lower left", fontsize=6, frameon=False)
fig.savefig(f"{T}/gene_tree_{aln}.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{T}/gene_tree_{aln}.svg", bbox_inches="tight")
Phylo.write(tree, f"{T}/gene_tree_{aln}.rooted.nwk", "newick")
print("tips drawn", len(leaves))
