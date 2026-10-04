#!/usr/bin/env python3
"""MAT locus size by idiomorph and genus, mapped onto the species tree.

Locus size = inner ends of the flanking genes (curator ruling 2026-10-03), only
for loci with polished flank genes on both sides of the core gene (calls.tsv
flank_status == both_sides). Held-out genomes (taxon_overrides 'hold out') are
dropped; confirmed overrides are counted under their likely genus.

Species tree: nf_phyling mucoromycota_odb12 FastTree tree (941 taxa), rooted on
Umbelopsis, pruned to one tip per genus that has >= MIN_LOCI sized loci. The tip
for a genus is a reference genome (BFD/Jena) when one exists, else an LCG genome
not flagged by the B12 name check.

Outputs: locus_size_by_genus.tsv (n, median, IQR per genus x idiomorph, with the
flank pairs seen), locus_size_loci.tsv (one row per sized locus), and
locus_size_tree.{png,svg}.
"""
import csv, re, statistics
from collections import defaultdict, Counter
from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat"
O = f"{C}/locus_size"
WT = "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-7c7ed99"
SPTREE = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-01_lcg_name_check/nf_out/protein/buildtree/mucoromycota_odb12/fasttree/protein-lcg_namecheck_v1-taxa_941.mucoromycota_odb12.fasttree.support.treefile"
FLAGS = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-01_lcg_name_check/flags.tsv"
TAXA = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-01_lcg_name_check/taxa.tsv"
MIN_LOCI = 2
COL = {"Plus": "#d62728", "Minus": "#1f77b4"}

ov = {}
for l in open(f"{WT}/db/taxon_overrides.tsv"):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 9 and not l.startswith("#") and f[0] != "genome_id":
        ov[f[0]] = f

def genus_of(r):
    o = ov.get(r["genome"])
    if o and o[5] == "confirmed":
        return o[2].split()[0]
    return r["name"].split()[0]

loci = []
for r in csv.DictReader(open(f"{C}/calls.tsv"), delimiter="\t"):
    if r["status"] != "called" or r["flank_status"] != "both_sides" or r["idiomorph"] not in COL:
        continue
    o = ov.get(r["genome"])
    if o and ("hold out" in o[6] or o[5] == "unconfirmed"):
        continue
    r["genus"] = genus_of(r)
    r["size"] = int(r["locus_size_inner"])
    r["pair"] = "|".join(sorted([r["left_flank"], r["right_flank"]]))
    loci.append(r)
with open(f"{O}/locus_size_loci.tsv", "w") as fo:
    keys = ["genus", "name", "source", "genome", "idiomorph", "confidence", "locus_class", "contig",
            "locus_start", "locus_end", "left_flank", "right_flank", "pair", "size"]
    w = csv.DictWriter(fo, fieldnames=keys, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    w.writeheader(); w.writerows(sorted(loci, key=lambda r: (r["genus"], r["idiomorph"], r["size"])))

by = defaultdict(list)
for r in loci:
    by[(r["genus"], r["idiomorph"])].append(r)
genera = sorted({g for g, _ in by if sum(len(by.get((g, i), [])) for i in COL) >= MIN_LOCI})

def q(v, p):
    v = sorted(v); k = (len(v) - 1) * p; lo = int(k); hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (k - lo)

rows = []
for g in genera:
    row = {"genus": g}
    for i in COL:
        v = [r["size"] for r in by.get((g, i), [])]
        row[f"{i}_n"] = len(v)
        row[f"{i}_median"] = int(statistics.median(v)) if v else ""
        row[f"{i}_q1"] = int(q(v, .25)) if v else ""
        row[f"{i}_q3"] = int(q(v, .75)) if v else ""
        row[f"{i}_min"] = min(v) if v else ""
        row[f"{i}_max"] = max(v) if v else ""
    pairs = Counter(r["pair"] for i in COL for r in by.get((g, i), []))
    row["flank_pairs"] = "; ".join(f"{k} ({n})" for k, n in pairs.most_common())
    rows.append(row)

# species tree: one tip per genus
taxa = {r["label"]: r for r in csv.DictReader(open(TAXA), delimiter="\t")}
tree = Phylo.read(SPTREE, "newick")
umb = [t for t in tree.get_terminals() if "Umbelopsis" in t.name]
tree.root_with_outgroup(tree.common_ancestor(umb))
lcg_name = {r["genome"]: genus_of(r) for r in csv.DictReader(open(f"{C}/calls.tsv"), delimiter="\t") if r["source"] == "LCG"}
# One tip per genus. Prefer reference tips (BFD, Jena), then LCG tips the B12 name
# check did not flag, so a misnamed genome never stands for its genus.
flagged = {r["genome"] for r in csv.DictReader(open(FLAGS), delimiter="\t")}
cands = {}
for t in tree.get_terminals():
    lab = t.name
    if lab.startswith("LCG__"):
        gid = lab[5:]
        if gid in flagged:
            continue
        g, rank = (lcg_name.get(gid) or gid.split("_")[0]), 1
    else:
        tr = taxa.get(lab, {})
        g = (tr.get("curator_taxonomy") or tr.get("genus") or re.sub(r"^(JENA|BFD)__", "", lab)).split()[0].split("_")[0]
        rank = 0
    if g in genera and (g not in cands or rank < cands[g][0]):
        cands[g] = (rank, t.name)
rep = {g: v[1] for g, v in cands.items()}
missing = [g for g in genera if g not in rep]
for t in list(tree.get_terminals()):
    if t.name not in rep.values():
        tree.prune(t)
tree.ladderize(reverse=True)
inv = {v: k for k, v in rep.items()}
order = [inv[t.name] for t in tree.get_terminals()]

with open(f"{O}/locus_size_by_genus.tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    for g in order + missing:
        w.writerow(next(r for r in rows if r["genus"] == g))

# figure
n = len(order)
fig, (axt, axd) = plt.subplots(1, 2, figsize=(11, 0.32 * n + 1.5), sharey=True,
                               gridspec_kw={"width_ratios": [1, 2.4], "wspace": 0.02})
depth = tree.depths(unit_branch_lengths=False)
ypos = {t: i for i, t in enumerate(tree.get_terminals())}
def y_of(c):
    if c.is_terminal():
        return ypos[c]
    ys = [y_of(k) for k in c.clades]
    c._y = (min(ys) + max(ys)) / 2
    return c._y
y_of(tree.root)
def draw(c):
    x, y = depth[c], (ypos[c] if c.is_terminal() else c._y)
    for k in c.clades:
        xk, yk = depth[k], (ypos[k] if k.is_terminal() else k._y)
        axt.plot([x, x], [y, yk], color="black", lw=0.8)
        axt.plot([x, xk], [yk, yk], color="black", lw=0.8)
        draw(k)
draw(tree.root)
xmax = max(depth.values())
for t, i in ypos.items():
    axt.text(depth[t] + xmax * 0.02, i, inv[t.name], va="center", fontsize=8, style="italic")
axt.set_xlim(0, xmax * 1.6); axt.axis("off")
for i, g in enumerate(order):
    for k, idio in enumerate(COL):
        v = [r["size"] / 1000 for r in by.get((g, idio), [])]
        if not v:
            continue
        dy = -0.17 if idio == "Plus" else 0.17
        axd.scatter(v, [i + dy] * len(v), s=10, color=COL[idio], alpha=0.6, lw=0)
        axd.plot([statistics.median(v)] * 2, [i + dy - 0.15, i + dy + 0.15], color="black", lw=1.2)
    nP, nM = len(by.get((g, "Plus"), [])), len(by.get((g, "Minus"), [])),
    axd.text(1.02, i, f"{nP}/{nM}", transform=axd.get_yaxis_transform(), va="center", fontsize=7)
axd.set_xscale("log")
axd.set_xlabel("MAT locus size, kb (inner ends of flanking genes; log scale)")
axd.text(1.02, -1.0, "n P/M", transform=axd.get_yaxis_transform(), fontsize=7)
axd.set_ylim(n - 0.5, -0.8)
axd.grid(axis="x", ls=":", lw=0.5)
axd.scatter([], [], color=COL["Plus"], label="Plus (sexP)"); axd.scatter([], [], color=COL["Minus"], label="Minus (sexM)")
axd.legend(loc="lower right", fontsize=8, frameon=False)
fig.savefig(f"{O}/locus_size_tree.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{O}/locus_size_tree.svg", bbox_inches="tight")
print("loci", len(loci), "genera", len(genera), "on tree", len(order), "missing from tree", missing)
