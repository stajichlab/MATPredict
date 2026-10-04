#!/usr/bin/env python3
"""Where do the non-MAT HMG genes (gene-tree outgroup) sit, and are they missed sexP/sexM?

Outgroup = classifier paralog negatives (non-locus HMG copies that fell outside the
sexP/sexM clades of the 2026-09-27 FastTree HMG-box tree) + P1. For each outgroup tip
in the 2026-10-03 gene trees (full length, HMG box):
- smallest clade holding the tip and >= 1 MAT tip: MAT types in it (P, M, both),
  UFBoot of that node, its size;
- nearest MAT tip by patristic distance and its type;
- genome context from the 2026-09-27 candidate table (contig, start, end, species) and
  the v0.6.0 calls (calls.tsv): does the genome have a Plus/Minus call, same contig,
  distance to it; the 2026-09-27 classifier verdict and margin.
Writes hmg_outgroup_placement.tsv; prints summaries.
"""
import csv, collections
from Bio import Phylo
C = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat"
T = f"{C}/tree"
O = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-04_unclassified_hmg"
tips = {r["tip"]: r for r in csv.DictReader(open(f"{T}/tips.tsv"), delimiter="\t")}
cand = {r["candidate_id"]: r for r in csv.DictReader(open(
    "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_sexMP_fasttree/candidates.tsv"), delimiter="\t")}
calls = collections.defaultdict(list)
for r in csv.DictReader(open(f"{C}/calls.tsv"), delimiter="\t"):
    if r["status"] == "called":
        calls[r["genome"]].append(r)
# cd-hit clusters: an outgroup representative may stand for several proteins
def clusters(path):
    out, cur = {}, None
    for l in open(path):
        if l.startswith(">Cluster"):
            cur = []
        else:
            name = l.split(">", 1)[1].split("...", 1)[0]
            cur.append(name)
            if l.rstrip().endswith("*"):
                rep = name
            out.setdefault("__members__", []).append((name, cur))
    return out
def sup(c):
    n = str(c.name or "")
    return float(n.split("/")[-1]) if "/" in n else None
def kind(x):
    return tips[x.name]["idiomorph"]
rows = []
for aln in ("full", "hmg"):
    t = Phylo.read(f"{T}/gene_tree_{aln}.rooted.nwk", "newick")
    L = t.get_terminals()
    mat = [x for x in L if kind(x) in ("Plus", "Minus")]
    parents = {}
    for cl in t.find_clades(order="level"):
        for ch in cl.clades:
            parents[ch] = cl
    for x in L:
        if kind(x) != "outgroup":
            continue
        node = x
        while node in parents:
            node = parents[node]
            ks = collections.Counter(kind(y) for y in node.get_terminals())
            if ks["Plus"] or ks["Minus"]:
                break
        types = "both" if ks["Plus"] and ks["Minus"] else ("P" if ks["Plus"] else "M")
        near = min(mat, key=lambda y: t.distance(x, y))
        src = tips[x.name]["label"].split(" ")[0]
        cr = cand.get(src, {})
        g = cr.get("genome") or tips[x.name]["genome"]
        gc = calls.get(g, [])
        same = [c for c in gc if c["contig"] == cr.get("contig")]
        dist = ""
        if same and cr.get("start"):
            s, e = int(cr["start"]), int(cr["end"])
            dist = min(max(0, int(c["locus_start"]) - e, s - int(c["locus_end"])) for c in same)
        rows.append(dict(tree=aln, tip=x.name, source=src, genome=g, species=cr.get("species", ""),
                         family=cr.get("family", ""), contig=cr.get("contig", ""), start=cr.get("start", ""),
                         end=cr.get("end", ""), clade_types=types, clade_ufboot=sup(node),
                         clade_mat_tips=ks["Plus"] + ks["Minus"], clade_outgroup_tips=ks["outgroup"],
                         nearest_mat=kind(near), nearest_dist=round(t.distance(x, near), 3),
                         genome_calls="|".join(sorted({c["idiomorph"] for c in gc})) or "none",
                         same_contig_as_call=bool(same), dist_to_call_bp=dist,
                         clf_2709=cr.get("clf_verdict", ""), clf_margin_2709=cr.get("clf_margin", ""),
                         clade_2709=cr.get("tree_clade", "")))
with open(f"{O}/hmg_outgroup_placement.tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
Cn = collections.Counter
for aln in ("full", "hmg"):
    R = [r for r in rows if r["tree"] == aln]
    print(f"\n== {aln}: {len(R)} outgroup tips")
    print(" first clade with MAT tips:", Cn(r["clade_types"] for r in R))
    for ty in ("P", "M"):
        rr = [r for r in R if r["clade_types"] == ty]
        strong = [r for r in rr if (r["clade_ufboot"] or 0) >= 95]
        print(f"  joins a {ty}-only clade: {len(rr)}; with UFBoot>=95: {len(strong)}; "
              f"genome calls {Cn(r['genome_calls'] for r in rr).most_common(4)}; "
              f"same contig as a call {sum(r['same_contig_as_call'] for r in rr)}")
    print(" nearest MAT tip:", Cn(r["nearest_mat"] for r in R))
    print(" classifier verdict (2026-09-27):", Cn(r["clf_2709"] for r in R).most_common(5))
