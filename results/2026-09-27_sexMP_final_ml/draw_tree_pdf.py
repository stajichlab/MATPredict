"""Draw the HMG-box tree as a tall, readable PDF.

Tip label = species + strain (from BFD samples.csv; for curated references,
the record's species + strain). Tip colour = detection label (Plus / Minus /
none) or reference type. Filled node dots = UFBoot >= 95; open = 80-94.
Rooted on the Ascomycota MAT1-2-1 outgroup, ladderized.

Usage: python3 draw_tree_pdf.py [PREFIX]   (default tree_hmm) -> PREFIX.pdf
"""
import csv, glob, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import yaml
from Bio import Phylo

prefix = sys.argv[1] if len(sys.argv) > 1 else "iq"
STRIP_P = sys.argv[2] if len(sys.argv) > 2 else "sexP clade"
STRIP_M = sys.argv[3] if len(sys.argv) > 3 else "sexM group"
SUPLAB = sys.argv[4] if len(sys.argv) > 4 else "UFBoot"
TITLE = sys.argv[5] if len(sys.argv) > 5 else prefix
DB = "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts/db"
SAMPLES = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"

tips = {r["tip"]: r for r in csv.DictReader(open("tip_names.tsv"), delimiter="\t")}
samp = {r["ASMID"]: r for r in csv.DictReader(open(SAMPLES))}
rec = {}
for f in glob.glob(f"{DB}/*/*/*/metadata.yaml"):
    try:
        o = yaml.safe_load(open(f))["organism"]
        rec[f.split("/")[-2]] = (o.get("species", ""), str(o.get("strain", "") or ""))
    except Exception:
        pass

COL = {"REF_sexM": "#1f4fbf", "REF_sexP": "#c41e1e", "REF_MAT1-2-1": "#6a3d9a",
       "Minus": "#2b8cbe", "Plus": "#e6550d", "-": "#555555", "undetermined": "#555555"}


def tip_text(r):
    if r["status"] == "reference":
        src = r["source"].split("|")
        rid = src[2] if src[0] == "OUT" else src[1]
        sp, st = rec.get(rid, (r["species"], ""))
        return f"{sp} {st}  [{r['kind'].replace('REF_', 'REF ')}]".replace("  ", " ", 0)
    genome = r["source"].split("|")[0]
    s = samp.get(genome, {})
    sp = s.get("SPECIES") or s.get("SPECIES_IN") or r["species"] or genome
    st = s.get("STRAIN", "")
    extra = {"call": "", "withheld": " (withheld)", "nonlocus": " (non-locus copy)"}.get(r["status"], "")
    k = int(r["n_collapsed"] or 1)
    dup = f" x{k}" if k > 1 else ""
    return f"{sp} {st}{extra}{dup}".strip()


tree = Phylo.read(f"{prefix}.treefile", "newick")
import os
con = Phylo.read(f"{prefix}.contree" if os.path.exists(f"{prefix}.contree") else f"{prefix}.treefile", "newick")
# support from the consensus tree, matched by tip set
sup = {}
for cl in con.find_clades():
    if not cl.is_terminal() and cl.confidence is not None:
        c = cl.confidence
        sup[frozenset(t.name for t in cl.get_terminals())] = c * 100 if c <= 1 else c
# Root on the Fusarium graminearum MAT1-2-1 (record 5518_3639_MAT_combined).
# The two MAT1-2-1 outgroups are NOT monophyletic in this tree, so rooting on
# both is undefined; this one splits the tree into a sexP half and a sexM half.
ROOT_RECORD = "5518_3639_MAT_combined"
out = [t for t in tree.get_terminals()
       if tips[t.name]["kind"] == "REF_MAT1-2-1" and ROOT_RECORD in tips[t.name]["source"]]
assert len(out) == 1, out
tree.root_with_outgroup(out[0])
tree.ladderize(reverse=True)

terms = tree.get_terminals()
n = len(terms)
ypos = {t: i for i, t in enumerate(terms)}
depth = tree.depths()
if not max(depth.values()):
    depth = tree.depths(unit_branch_lengths=True)


def y_of(cl):
    if cl in ypos:
        return ypos[cl]
    ys = [y_of(c) for c in cl.clades]
    ypos[cl] = (min(ys) + max(ys)) / 2
    return ypos[cl]


y_of(tree.root)
maxx = max(depth.values())
row_in = 0.13
fig_h = max(8, n * row_in + 2)
fig, ax = plt.subplots(figsize=(14, fig_h))
lw = 0.5
for cl in tree.find_clades():
    x, y = depth[cl], ypos[cl]
    if cl.clades:
        ys = [ypos[c] for c in cl.clades]
        ax.plot([x, x], [min(ys), max(ys)], color="black", lw=lw)
        for c in cl.clades:
            ax.plot([x, depth[c]], [ypos[c], ypos[c]], color="black", lw=lw)
        s = sup.get(frozenset(t.name for t in cl.get_terminals()))
        if s is not None and cl != tree.root:
            if s >= 95:
                ax.plot(x, y, "o", ms=2.2, color="black")
            elif s >= 80:
                ax.plot(x, y, "o", ms=2.2, mfc="white", mec="black", mew=0.4)
for t in terms:
    r = tips[t.name]
    key = r["kind"] if r["status"] == "reference" else r["idiomorph_label"]
    c = COL.get(key, "#555555")
    bold = r["status"] == "reference"
    ax.text(depth[t] + maxx * 0.005, ypos[t], tip_text(r), va="center", ha="left",
            fontsize=5.2, color=c, fontweight="bold" if bold else "normal")
# clade strip at the right edge: tree_clade from read_tree.py (refs by their type)
strip_x = maxx * 1.42
CLADE_COL = {"sexM": "#1f4fbf", "sexP": "#c41e1e"}
for t in terms:
    r = tips[t.name]
    cl = r["tree_clade"] or {"REF_sexM": "sexM", "REF_sexP": "sexP"}.get(r["kind"], "")
    if cl in CLADE_COL:
        ax.plot([strip_x, strip_x], [ypos[t] - 0.5, ypos[t] + 0.5], color=CLADE_COL[cl], lw=4,
                solid_capstyle="butt")
ax.text(strip_x, -1.5, "clade", ha="center", fontsize=6)
handles_strip = [Line2D([], [], color=CLADE_COL["sexM"], lw=4, label=STRIP_M),
                 Line2D([], [], color=CLADE_COL["sexP"], lw=4, label=STRIP_P)]
ax.set_ylim(n, -1)
ax.set_xlim(0, maxx * 1.45)
ax.axis("off")
ax.plot([0, 0.5], [n - 0.2, n - 0.2], color="black", lw=0.8)
ax.text(0.25, n + 0.8, "0.5 substitutions/site", ha="center", fontsize=6)
handles = [Line2D([], [], color=COL[k], lw=0, marker="s", ms=6, label=l) for k, l in [
    ("REF_sexM", "curated sexM reference"), ("REF_sexP", "curated sexP reference"),
    ("REF_MAT1-2-1", "Ascomycota MAT1-2-1 (root: F. graminearum)"), ("Minus", "detection label Minus"),
    ("Plus", "detection label Plus"), ("-", "no label (non-locus copy / undetermined)")]]
handles += handles_strip + [Line2D([], [], color="black", lw=0, marker="o", ms=4, label=f"{SUPLAB} >= 95"),
            Line2D([], [], mfc="white", mec="black", lw=0, marker="o", ms=4, label=f"{SUPLAB} 80-94")]
ax.legend(handles=handles, loc="upper left", fontsize=7, frameon=False, bbox_to_anchor=(0, 1.0))
ax.set_title(TITLE,
             fontsize=9, loc="left")
fig.savefig(f"{prefix}.pdf", bbox_inches="tight")
print(f"wrote {prefix}.pdf: {n} tips, {fig_h:.0f} in tall")
