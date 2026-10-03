"""Write copies of the tree with readable tip names, plus a FigTree colour file.

Tip name = species | family | what the leaf is | detection label | best-hit ref
gene | tree clade | leaf id.  Reference tips name the curated record and its
gene (sexM / sexP / MAT1-2-1).  "xN" = N identical HMG domains collapsed onto
this leaf (dedup_map.tsv).

Usage: python3 name_tree.py [PREFIX]   (default tree_hmm)
Writes PREFIX.named.treefile, PREFIX.named.contree, PREFIX.named.nex (FigTree,
tips coloured), tip_names.tsv.
"""
import collections, csv, glob, re, sys
import yaml

prefix = sys.argv[1] if len(sys.argv) > 1 else "tree_ft"
DB = "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts/db"

ids = dict(l.rstrip("\n").split("\t") for l in list(open("taxon_ids.tsv"))[1:])
cand = {r["candidate_id"]: r for r in csv.DictReader(open("candidates.tsv"), delimiter="\t")}
nmem = collections.Counter(r["rep"] for r in csv.DictReader(open("dedup_map.tsv"), delimiter="\t"))

record_species = {}
for f in glob.glob(f"{DB}/*/*/*/metadata.yaml"):
    rec = f.split("/")[-2]
    try:
        record_species[rec] = yaml.safe_load(open(f))["organism"]["species"]
    except Exception:
        pass


def safe(s):
    return re.sub(r"[^A-Za-z0-9._+-]+", "_", s).strip("_")


COLOUR = {"REF_sexM": "#1f4fbf", "REF_sexP": "#c41e1e", "REF_MAT1-2-1": "#6a3d9a",
          "Minus": "#6baed6", "Plus": "#fb6a4a", "other": "#7f7f7f"}
names, colour, rows = {}, {}, []
for tid, n in ids.items():
    if n.startswith(("REF|", "OUT|")):
        parts = n.split("|")
        rec = parts[2] if n.startswith("OUT|") else parts[1]
        gene = "MAT1-2-1" if n.startswith("OUT|") else parts[-1]
        kind = "REF_MAT1-2-1" if n.startswith("OUT|") else f"REF_{gene}"
        sp = record_species.get(rec, rec)
        label = f"{kind}|{sp}|{rec}|{tid}"
        col = COLOUR.get(kind, COLOUR["other"])
        row = dict(tip=tid, kind=kind, species=sp, family="", status="reference",
                   idiomorph_label="", best_ref_gene=gene, tree_clade="", n_collapsed=1, source=n)
    else:
        r = cand[n]
        ref_gene = r["best_ref"].split("|")[-1] if r["best_ref"] else ""
        status = "call" if r["status"] == "called" else r["status"]
        lab = r["idiomorph_label"] or "-"
        k = nmem.get(n, 1)
        label = "|".join([r["species"] or r["genome"], r["family"] or r["group"] or "?",
                          status, lab, f"hit_{ref_gene}", r["tree_clade"],
                          f"x{k}" if k > 1 else "", tid])
        col = COLOUR.get(lab, COLOUR["other"])
        row = dict(tip=tid, kind=r["kind"], species=r["species"], family=r["family"], status=status,
                   idiomorph_label=lab, best_ref_gene=ref_gene, tree_clade=r["tree_clade"],
                   n_collapsed=k, source=n)
    names[tid] = safe(label.replace("||", "|"))
    colour[tid] = col
    rows.append(row)

import os
for suffix in ("treefile", "contree"):
    if not os.path.exists(f"{prefix}.{suffix}"):
        continue
    txt = open(f"{prefix}.{suffix}").read()
    txt = re.sub(r"\b(T\d{4})\b(?=[:,)])", lambda m: names[m.group(1)], txt)
    open(f"{prefix}.named.{suffix}", "w").write(txt)

con = open(f"{prefix}.named.contree" if os.path.exists(f"{prefix}.named.contree") else f"{prefix}.named.treefile").read().strip()
with open(f"{prefix}.named.nex", "w") as fo:
    fo.write(f"#NEXUS\nbegin taxa;\n\tdimensions ntax={len(names)};\n\ttaxlabels\n")
    for tid in sorted(names):
        fo.write(f"\t'{names[tid]}'[&!color={colour[tid]}]\n")
    fo.write(";\nend;\n\nbegin trees;\n\ttree tree_1 = [&U] " + con + "\nend;\n")

with open("tip_names.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=["name"] + list(rows[0]), delimiter="\t")
    w.writeheader()
    for row in rows:
        w.writerow({"name": names[row["tip"]], **row})
print(f"{len(names)} tips named; wrote {prefix}.named.treefile/.contree/.nex and tip_names.tsv")
print(collections.Counter(r["kind"] if r["status"] == "reference" else r["idiomorph_label"] for r in rows))
