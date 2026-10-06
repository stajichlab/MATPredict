#!/usr/bin/env python3
"""Merge the per-genome locus tables (build_genome.py) with labels and provenance.

Label rules
 - Agaricomycete panel genomes (6 of PR #33 + 4 new): label from panel_loci.tsv by overlap
   with the curated B-locus receptor cluster (same rule as PR #33).
 - All other genomes: a locus is a mating receptor if its best curated-reference hit is the
   genome's OWN db record receptor with identity >= 0.95 (the record protein is the genome's
   own gene, so identity is ~1.0). Every other STE3-like locus is "other".
Writes labelled_loci.tsv and labelled.faa (protein per locus, key = locus_id).
Usage: make_labels.py OUTDIR PANEL_LOCI.tsv
"""
import csv
import glob
import os
import sys

OUT, PANEL = sys.argv[1], sys.argv[2]

# asm prefix -> (species, order, own record id prefix or "", grade of the mating label, evidence, headline)
G = {
    "GCA_016772295.1": ("Coprinopsis cinerea", "Agaricales", "", "B", "B43 locus (PMID 9539426, 10757757); mapped, per-gene specificity not re-verified here", True),
    "GCF_000143185.2": ("Schizophyllum commune", "Agaricales", "", "A", "bar3/bbr2 B-locus receptors, transformation tests (PMID 7489716)", True),
    "GCF_000271585.1": ("Trametes versicolor", "Polyporales", "", "C", "genome-derived record, chosen from CAAX positional evidence (tier 2)", False),
    "GCA_001683735.1": ("Grifola frondosa", "Polyporales", "", "C", "genome-derived record, chosen from CAAX positional evidence (tier 2)", False),
    "GCA_984573805.1": ("Russula nobilis", "Russulales", "", "C", "genome-derived record, chosen from CAAX positional evidence (tier 2)", False),
    "GCF_000320585.1": ("Heterobasidion irregulare", "Russulales", "", "C", "genome-derived record, chosen from CAAX positional evidence (tier 2)", False),
    "GCF_000328475.2": ("Ustilago maydis 521", "Ustilaginales", "", "A", "pra1 at the a locus (PMID 1310895, Bolker 1992; U37795)", True),
    "GCF_000091045.1": ("Cryptococcus neoformans JEC21", "Tremellales", "", "A", "STE3alpha at MAT alpha (PMID 12455690)", True),
    "GCA_056621545.1": ("Cryptococcus neoformans JEC20", "Tremellales", "40410_jec20_MAT_a", "A", "STE3a at MAT a (PMID 12455690)", True),
    "GCA_000988875.2": ("Rhodotorula toruloides NBRC0880", "Sporidiobolales", "", "B", "STE3a2 at the P/R locus (biorxiv 10.1101/2025.09.11.675505)", True),
    "GCA_921037615.3": ("Rhodotorula toruloides CBS14", "Sporidiobolales", "", "B", "STE3a1 at the P/R locus (biorxiv 10.1101/2025.09.11.675505)", True),
    "GCA_026119225.1": ("Rhodotorula toruloides JJ10-1", "Sporidiobolales", "29898_jj10-1_redPR_A2", "B", "STE3a2 at the P/R locus", True),
    "GCA_920103745.3": ("Rhodotorula toruloides CBS20", "Sporidiobolales", "5535_cbs-20_redPR_A1", "B", "STE3a1 at the P/R locus", True),
    "GCA_024748845.1": ("Rhodotorula toruloides JY1105", "Sporidiobolales", "5537_jy1105_redPR_A1", "B", "STE3a1 at the P/R locus", True),
    "GCA_023212685.1": ("Ustilaginales CBS10937", "Ustilaginales", "1652704_cbs-10937_aLocus_a1", "B", "pra1 at the a locus (PMID 37312063)", True),
    "GCA_023212835.1": ("Ustilaginales CBS10006", "Ustilaginales", "203535_cbs-10006_aLocus_a1", "B", "pra1 at the a locus (PMID 37312063)", True),
    "GCA_023212605.1": ("Ustilaginales CBS10005", "Ustilaginales", "203536_cbs-10005_aLocus_a1", "B", "pra1 at the a locus (PMID 37312063)", True),
    "GCA_023212725.1": ("Ustilaginales CBS131475", "Ustilaginales", "349360_cbs-131475_aLocus_a1", "B", "pra1 at the a locus (PMID 37312063)", True),
    "GCA_023212615.1": ("Ustilaginales CBS131463", "Ustilaginales", "49012_cbs-131463_aLocus_a1", "B", "pra1 at the a locus (PMID 37312063)", True),
    "GCA_023212695.1": ("Ustilaginales CBS10417", "Ustilaginales", "63387_cbs-10417_aLocus_a1", "B", "pra1 at the a locus (PMID 37312063)", True),
    "GCA_023212635.2": ("Ustilaginales CBS167882", "Ustilaginales", "84751_cbs-167882_aLocus_a1", "B", "pra1 at the a locus (PMID 37312063)", True),
    "GCA_056320075.1": ("Wallemia EXF-10342", "Wallemiales", "1708542_exf-10342_wallMAT_v2", "C", "PUTATIVE Wallemia MAT locus (PMID 22326418, 31167502)", False),
    "GCF_000263375.1": ("Wallemia mellicola CBS633.66", "Wallemiales", "671144_cbs-633-66_wallMAT_v1", "C", "PUTATIVE Wallemia MAT locus (PMID 22326418, 31167502)", False),
}
PANEL_ASM = set()
panel = []
for r in csv.DictReader(open(PANEL), delimiter="\t"):
    panel.append(r)
    PANEL_ASM.add(r["asm"])
rows, fa = [], open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "labelled.faa"), "w")
prot = {}
for f in glob.glob(os.path.join(OUT, "*.prot.faa")):
    k = None
    for line in open(f):
        if line.startswith(">"):
            k = line[1:].strip()
            prot[k] = ""
        else:
            prot[k] += line.strip()
XO = {}
for f in glob.glob(os.path.join(OUT, "*.xo.tsv")):
    for x in csv.DictReader(open(f), delimiter="\t"):
        XO[(x["asm"], x["contig"], x["start"], x["end"])] = x
for f in sorted(glob.glob(os.path.join(OUT, "*.loci.tsv"))):
    for r in csv.DictReader(open(f), delimiter="\t"):
        x = XO[(r["asm"], r["contig"], r["start"], r["end"])]
        for k in ("d_Hx_xo", "nHx_10kb_xo", "d_HD_xo", "d_FLANK_xo", "n_HD_20kb_xo", "n_FLANK_50kb_xo"):
            r[k] = x[k]
        asm = r["asm"]
        pre = asm[:15]
        sp, order, own, grade, ev, headline = G[pre]
        s, e = int(r["start"]), int(r["end"])
        lab = None
        if asm in PANEL_ASM:
            hit = [p for p in panel if p["asm"] == asm and p["contig"] == r["contig"] and int(p["start"]) <= e and int(p["end"]) >= s]
            if hit:
                lab = hit[0]["mating"]
        else:
            ref = r["best_ref"]
            own_hit = own and ref.startswith("REF|" + own) and float(r["best_ref_ident"]) >= 0.95
            lab = "mating" if own_hit else "other"
        if lab is None:
            lab = "unlabelled"
        lid = f"{asm}|{r['contig']}:{s}-{e}"
        r["locus_id"] = lid
        r["species"], r["order"], r["label"] = sp, order, lab
        r["grade"] = grade if lab == "mating" else ("D" if lab == "other" else "NA")
        r["evidence"] = ev
        r["headline"] = headline
        r["has_prot"] = int(lid in prot)
        rows.append(r)
        if lid in prot:
            fa.write(f">{lid}\n{prot[lid]}\n")
cols = ["locus_id", "species", "order", "label", "grade", "headline", "evidence"] + [c for c in rows[0] if c not in ("locus_id", "species", "order", "label", "grade", "headline", "evidence")]
with open(os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "labelled_loci.tsv"), "w") as fo:
    fo.write("\t".join(cols) + "\n")
    for r in rows:
        fo.write("\t".join(str(r[c]) for c in cols) + "\n")
from collections import Counter
print(Counter((r["order"], r["label"]) for r in rows))
print(Counter(r["has_prot"] for r in rows))
