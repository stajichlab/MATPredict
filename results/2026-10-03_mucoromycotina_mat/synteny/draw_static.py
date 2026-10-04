#!/usr/bin/env python3
"""Static synteny figure (pyGenomeViz) of the 26 selected loci.

Tracks: one per locus, in species-tree order (selection.tsv), each oriented so
its core MAT gene points right (reverse-complemented when needed). Genes coloured by
MATPredict name (sexP, sexM, flanks; other genes grey). Links between adjacent
tracks: clinker's gene-to-gene protein identities (clinker_alignments.txt),
shaded by identity (>= 30%). Each track is centred on its core MAT gene so the
idiomorphs line up.
Writes mucoromycotina_MAT_synteny.{png,svg}.
"""
import csv, re
from Bio import SeqIO
from pygenomeviz import GenomeViz
from matplotlib.patches import Patch

O = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat/synteny"
COL = {"sexP": "#d62728", "sexM": "#1f77b4", "tptA": "#2ca02c", "rnhA": "#9467bd", "glrA": "#ff7f0e",
       "algA": "#8c564b", "btbA": "#e377c2", "other": "#c7c7c7"}
sel = list(csv.DictReader(open(f"{O}/selection.tsv"), delimiter="\t"))
recs, genes = {}, {}
for s in sel:
    r = SeqIO.read(f"{O}/gbk/{s['sample']}.gbk", "genbank")
    core = "sexP" if s["idiomorph"] == "Plus" else "sexM"
    cf = [f for f in r.features if f.type == "CDS" and f.qualifiers["gene"][0] == core]
    if cf and cf[0].location.strand == -1:
        # orient every locus so its core MAT gene is on the forward strand
        r = r.reverse_complement(id=r.id, name=r.name, description=r.description,
                                 features=True, annotations=True)
        s["flipped"] = "yes"
    recs[s["sample"]] = r
    for f in r.features:
        if f.type == "CDS":
            genes[f.qualifiers["locus_tag"][0]] = (s["sample"], int(f.location.start), int(f.location.end),
                                                   f.location.strand, f.qualifiers["gene"][0])
# centre each track on its core gene
offsets, width = {}, 0
for s in sel:
    core = "sexP" if s["idiomorph"] == "Plus" else "sexM"
    cs = [g for g in genes.values() if g[0] == s["sample"] and g[4] == core]
    mid = (cs[0][1] + cs[0][2]) // 2 if cs else len(recs[s["sample"]].seq) // 2
    offsets[s["sample"]] = mid
    width = max(width, mid, len(recs[s["sample"]].seq) - mid)

gv = GenomeViz(fig_width=14, fig_track_height=0.45, track_align_type="left", feature_track_ratio=0.35)
gv.set_scale_bar(ymargin=0.5)
names = {}
for s in sel:
    label = f"{s['name']} ({s['idiomorph']}{', ' + s['source'] if s['source'] != 'LCG' else ''})"
    label = re.sub(r"\s+", " ", label)
    L = len(recs[s["sample"]].seq)
    off = width - offsets[s["sample"]]
    tr = gv.add_feature_track(label, segments=(0, L), offset=off, labelsize=9)
    names[s["sample"]] = label
    for f in recs[s["sample"]].features:
        if f.type != "CDS":
            continue
        g = f.qualifiers["gene"][0]
        tr.add_features(f, fc=COL.get(g, COL["other"]), ec="black", lw=0.3, plotstyle="arrow",
                        label_type=None)
# links from clinker
pair = None
for line in open(f"{O}/clinker_alignments.txt"):
    m = re.match(r"^(\S+) vs (\S+)$", line.strip())
    if m:
        pair = (m.group(1), m.group(2)); continue
    f = line.split()
    if len(f) == 4 and f[0] in genes and f[1] in genes:
        a, b = genes[f[0]], genes[f[1]]
        ia = [x["sample"] for x in sel].index(a[0]); ib = [x["sample"] for x in sel].index(b[0])
        if abs(ia - ib) != 1:
            continue
        ident = float(f[2]) * 100
        if ident < 30:
            continue
        gv.add_link((names[a[0]], a[1], a[2]), (names[b[0]], b[1], b[2]), color="grey",
                    inverted_color="tomato", v=ident, vmin=30, vmax=100, alpha=0.6)
fig = gv.plotfig()
gv.set_colorbar(["grey", "tomato"], vmin=30, vmax=100, bar_label="% identity", bar_labelsize=8)
fig.legend(handles=[Patch(fc=v, ec="black", lw=0.3, label=k) for k, v in COL.items()],
           loc="upper center", ncol=len(COL), fontsize=8, frameon=False, bbox_to_anchor=(0.5, 1.02))
fig.savefig(f"{O}/mucoromycotina_MAT_synteny.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{O}/mucoromycotina_MAT_synteny.svg", bbox_inches="tight")
print("ok", len(sel), "tracks")
