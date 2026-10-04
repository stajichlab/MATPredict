# sexP and sexM gene trees across Mucoromycotina (2026-10-04)

Status: decided (descriptive); RAxML-NG check of the HMG-box tree running.
Part of the Mucoromycotina MAT campaign
([2026-10-03 report](2026-10-03_mucoromycotina-mat-campaign.md)).

## Question
In a rooted tree of the sexP (Plus) and sexM (Minus) HMG proteins called in
BFD, LCG and Jena genomes:
1. Do sexP and sexM each form one clade, with no tip of the other type?
2. How well supported are they, and does the answer depend on the alignment
   (HMG box only versus full length)?
3. What other patterns are visible (genus grouping, Umbelopsis, divergence,
   non-MAT HMG genes inside the MAT clades, misnamed genomes)?

## Ruling
J. Stajich, 2026-10-03: root with a non-MAT HMG-box outgroup; label tips sexP
or sexM; build from LCG, Jena and BFD.

## Data
- Calls: v0.6.0 (`run-7c7ed99`), 973 genomes, 911 loci
  (`results/2026-10-03_mucoromycotina_mat/calls.tsv`).
- Ingroup (`tree/select_tree_set.py`): polished core proteins of called loci
  plus the curated record proteins. Held-out (`hold out`) genomes and proteins
  under 60 aa dropped. sexP 487, sexM 409 proteins.
- Outgroup: 190 non-MAT HMG proteins = the classifier's paralog-negative set
  (HMG copies away from the MAT locus that fell outside the sexP and sexM clades
  in the 2026-09-26/27 trees) plus the P1 sexM-like paralog.
- Redundancy: cd-hit 0.98 within each idiomorph (sexP 163, sexM 173
  representatives), 0.90 on the outgroup (190).
- Tip labels are the classifier calls (Plus/Minus), not curator labels.

## Method
| Step | HMG box | Full length |
|---|---|---|
| Alignment | hmmalign to Pfam PF00505, match columns only | MAFFT L-INS-i, ClipKIT kpic-smart-gap |
| Size | 523 sequences x 69 sites (68 parsimony-informative) | 526 x 526 (516 informative) |
| Model (ModelFinder, BIC; LG/WAG/JTT/VT/Q.pfam) | LG+R6 | JTT+F+R8 |
| Tree | IQ-TREE 3.0.1, UFBoot 1000, SH-aLRT 1000, seed 20261003 | same |
| Run time (16 CPUs) | 12 h 30 min (job 29389343) | 2 h 40 min (job 29380482) |

Rooting (`tree/draw_gene_tree.py`): on the outgroup tip farthest from the
ingroup, then on the outgroup clade if the outgroup is one clade. It is not one
clade in either tree, so both trees are rooted on one outgroup tip.

Clade test (`tree/tree_patterns.py`, output `tree/patterns.txt`): the "main
clade" of an idiomorph is the largest clade that holds no tip of the other
idiomorph (outgroup tips allowed).

## Results

### 1. Are sexP and sexM separate clades?
| | sexP | sexM |
|---|---|---|
| Full length | all 163 tips in one clade with no sexM tip, but that clade also holds 131 outgroup tips (UFBoot 29); MRCA of sexP tips alone: 165 tips (2 outgroup), UFBoot 40 | one clade, 173/173 tips + 1 outgroup tip, **UFBoot 93** |
| HMG box | one clade, 162/162 tips + 4 outgroup tips, **UFBoot 99** (sexP MRCA alone: UFBoot 93) | not one clade: largest sexM-only clade 48/173 (UFBoot 53); the sexP clade nests among the other sexM lineages |

- Where a clade is supported, it is clean: the full-length sexM clade holds no
  sexP tip, and the HMG-box sexP clade holds no sexM tip. No single sexP tip is
  placed among sexM tips, or the reverse, in either tree; the disagreement is
  only about the deep, weakly supported nodes.
- The two alignments give opposite answers on monophyly: the HMG box supports
  sexP and splits sexM; the full length supports sexM and, at its base, only
  weakly sexP. In the full-length tree, 152 of 163 sexP tips form a clade with
  no outgroup tip; the other 11 (Circinella 6, Rhizomucor 2, Thamnostylum,
  Phascolomyces, Absidia 1) sit among outgroup tips.
- The 69-site HMG box gives weak deep nodes. The sexM split there is
  unresolved; it is not evidence that sexP arose inside sexM.
- This repeats the 2026-09-26 trees on the same domain
  ([sexM/sexP trees](2026-09-26_sexMP-trees.md)): sexP clade UFBoot 99, sexM
  references not monophyletic (UFBoot 57). The new trees have about five times
  the tips and add the full-length alignment. Among this project's trees it is
  the first with all sexM in one clade at UFBoot >= 90 (the 2026-09-27 trimmed
  tree had the 9 sexM references together at UFBoot 74).

### 2. Rooting and the outgroup
- The non-MAT HMG outgroup is not monophyletic in either tree. In the full-
  length tree 132 outgroup tips fall inside the MRCA of all MAT tips (100 in the
  HMG-box tree), so sexP and sexM are not shown as sister groups, and the root
  position is uncertain.
- The outgroup set is not independent: it was defined as HMG copies outside the
  sexP/sexM clades of an earlier tree on the same domain.

### 3. Outgroup (non-MAT) HMG genes nested in the MAT clades
| Tip | Genome | In |
|---|---|---|
| GCA_025716815.1 Dicele1 h5 | Dichotomocladium elegans (Mucorales) | sexP (both trees) |
| GCA_016758965.1 h6 | Circinella minor (Mucorales) | sexP (both trees) |
| GCA_002261195.1 h3 | Bifiguratus adelaidae (Endogonales) | sexP (HMG box) |
| GCA_027478255.1 h4 | Dispira simplex (Dimargaritales) | sexP (HMG box) |
| GCF_025528875.1 Mycafr1 h1 | Mycotypha africana (Mucorales) | sexM (both trees) |

These were labelled non-MAT because they are not at a called locus. Three are
in Mucorales genomes; two are outside Mucorales (Endogonales, and Dimargaritales
in Kickxellomycotina). They may be MAT-gene copies away from the locus, or MAT
genes in lineages with no curated record. Not examined.

### 4. Genus grouping inside each idiomorph
| | sexP genera monophyletic | sexM genera monophyletic |
|---|---|---|
| Full length | 11/18 | 9/21 |
| HMG box | 10/18 | 10/21 |

Genera not monophyletic in both trees include Mucor (59 sexP / 53 sexM tips),
Rhizopus, Backusella, Absidia and Circinella. Genus names are curator names for
LCG/Jena and BFD names otherwise, so misnamed genomes add to this; the trees
were not compared with the species tree here.

### 5. Umbelopsis
- sexP: the 7 Umbelopsis tips form one clade in both trees. In the full-length
  tree it is sister to 145 of the other sexP tips (the remaining 11 are among
  outgroup tips), close to the species-tree position of Umbelopsidales as sister
  to Mucorales. In the HMG-box tree it is sister to Syncephalastrum (4 tips).
- sexM: the 6 Umbelopsis tips form one clade with one LCG genome named
  "Mucor sp. NRRL 1454", in both trees. Its sister is Backusella (full length)
  or Cunninghamella + Chaetocladium (HMG box): not the species-tree position.

### 6. A misnamed genome
LCG "Mucor sp. NRRL 1454" groups with Umbelopsis sexM in both gene trees. The
B12 name check placed it with Umbelopsis in both BUSCO species trees (support
1.0; `results/2026-10-01_lcg_name_check/flags.tsv`), but it has no row in
`db/taxon_overrides.tsv`. Curator, 2026-10-04: override added (Umbelopsis sp.,
genus, confirmed; LSU D1/D2 99.5% U. tibetica type; ANNOTATION_ERRORS C13).

### 7. Divergence
Median patristic distance between random tip pairs inside the main clade
(substitutions per site): full length sexP 3.76 (IQR 2.56-4.74), sexM 4.69
(3.58-5.50); HMG box sexP 1.09, sexM 0.99 (sexM main clade there is 48 tips).
On full length, sexM tips are further apart than sexP tips. This was not
tested for significance, and the sexP main clade there includes outgroup tips.

## Decision
Descriptive; no tool change. For figures: show the full-length tree as the
main tree (sexM one clade, UFBoot 93) and the HMG-box tree as a supplement, and
state that the root and the sexP-sexM relationship are not resolved.

## Open
- RAxML-NG check of the HMG-box tree (job 29397060).
- The 5 nested outgroup HMG genes: locus, gene order, identity.
- Compare genus placement with the species tree; test the sexM/sexP divergence
  difference.
- A better outgroup: HMG genes chosen by a criterion independent of the
  2026-09 sexP/sexM trees.

## Files
- `results/2026-10-03_mucoromycotina_mat/tree/`: inputs (`ingroup_sexP.faa`,
  `ingroup_sexM.faa`, `outgroup.faa`, `tips.tsv`, `*.clstr`), alignments
  (`hmg.afa`, `full.afa`), IQ-TREE outputs (`iq_hmg.*`, `iq_full.*`), rooted
  trees and figures (`gene_tree_{hmg,full}.{rooted.nwk,png,svg}`), summaries
  (`gene_tree_summary_{hmg,full}.tsv`, `patterns.txt`), scripts
  (`select_tree_set.py`, `run_tree.sh`, `run_iq_one.sh`, `draw_gene_tree.py`,
  `tree_patterns.py`).
- Figures: [full length](../results/2026-10-03_mucoromycotina_mat/tree/gene_tree_full.png),
  [HMG box](../results/2026-10-03_mucoromycotina_mat/tree/gene_tree_hmg.png).
