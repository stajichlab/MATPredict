# Non-MAT HMG genes inside the sexP/sexM gene trees: missed MAT genes?

Status: open (two candidates for curation; outgroup and negative set to revise).
Follows [sexP and sexM gene trees](2026-10-04_sexP-sexM-gene-trees.md).

## Question
Many genes of the non-MAT HMG outgroup fall within the sexP or sexM part of
the 2026-10-03 gene trees. Are they sexP or sexM genes that MATPredict did not
call, and can they be classified?

## Data
- Trees: `results/2026-10-03_mucoromycotina_mat/tree/gene_tree_{full,hmg}.rooted.nwk`.
- The outgroup is the classifier's paralog-negative set (189 proteins) plus P1:
  HMG copies away from any called locus that fell outside the sexP/sexM clades
  of the 2026-09-27 FastTree HMG-box tree. It is also the set that sets the
  MAT-gene gate threshold (95th percentile of their classifier scores). Most
  genomes contribute 2-3 copies; it includes Kickxellomycotina,
  Mortierellomycotina and Endogonales genomes.
- Genome calls: v0.6.0 campaign (`calls.tsv`).
- Script and tables: `results/2026-10-04_unclassified_hmg/`.

## Results

### 1. Most of the "many" is the rooting, not nested genes
The outgroup is not one clade, so most outgroup genes attach along the stem of
a MAT clade rather than inside it.
- Full length: 131 of 190 outgroup genes first meet MAT tips at a node whose
  MAT tips are all sexP, but that node sits below the sexP clade; only 2 are
  inside the sexP MRCA and 1 inside the sexM MRCA.
- HMG box: 83 of 188 first meet a sexM-only node, but only 3 such nodes have
  UFBoot >= 95; the sexM part of this tree is unresolved (69 sites).

### 2. A sexP-related HMG lineage (full-length tree)
The sexP clade (UFBoot 40) has a sister group of 40 non-MAT HMG genes; sexP plus
this group: UFBoot 100.
- Genomes: Mucorales 27 genes, Kickxellales 4, Dimargaritales 3, Mortierellales
  2, Umbelopsidales 2, Ramicandelaberales 1, Endogonales 1.
- Their genomes' MAT calls: Plus 14, Minus 9, both 1, none 16.
- Present in Plus and in Minus genomes, and next to a called sexP in Plus
  genomes: a conserved sexP-related paralog family, not a missed idiomorph gene.
  It is also the reason these genes sit close to sexP and score near the gate.
- Further down the sexP stem: a 72-gene group (mostly Mucorales, both mating
  types) and a 14-gene Mortierellales group (UFBoot 29 each).

### 3. Genes inside a supported MAT clade
| Gene | Genome | Position | Tree placement | Genome call | Flank genes within 60 kb |
|---|---|---|---|---|---|
| Mycafr1 h1 | Mycotypha africana NRRL 2978 (GCF_025528875.1) | NW_026515730.1:1,203,502-1,204,971 | sexM clade (both trees) | Plus (sexP at NW_026515730.1:1,357,549-1,365,893) | none |
| Dicele1 h5 | Dichotomocladium elegans RSA 919 (GCA_025716815.1) | JAIXMN010000002.1:75,748-76,284 (contig 2.42 Mb) | sexP, with the Circinella/Thamnostylum/Phascolomyces sexP (UFBoot 96 full, 95 HMG box) | none (uncalled) | none |
| Bifiguratus h3 | Bifiguratus adelaidae (Endogonales) | not located | sexP clade, HMG-box tree only (UFBoot 99) | not in the campaign | not checked |
| C. minor h6 | Circinella minor (GCA_016758965.1) | not located | inside the sexP MRCA (full length) | Plus | not checked |

- **Mycotypha africana is a published homothallic species with sexM about
  150 kb from sexP on the same scaffold (Schulz 2016).** h1 is on the scaffold
  of the called sexP, 153 kb from it, and groups with sexM in both trees. It is
  very likely that sexM. MATPredict does not call it, and the gene is in the
  classifier's negative set. No flank gene lies within 60 kb, so the rule that
  a locus needs two modelled genes probably withholds it (not checked in the
  report).
  tblastn also finds two more HMG hits next to it (83-89% identity, within 1.5 kb),
  not examined.
- Dichotomocladium elegans has no call. h5 groups with the Lichtheimiaceae-type
  sexP, away from any flank gene. Candidate missed sexP (or a sexP copy outside
  the locus); needs synteny and annotation review.

## Effects on MATPredict
- The negative set holds at least one probable MAT gene (Mycafr1 h1). The gate
  threshold is the 95th percentile of 189 negatives, so one gene changes it
  little, but it should be removed and the gate recomputed (classifier rebuild
  needs curator approval).
- Homothallic loci without flanks (Mycotypha) are not called by design. A
  report-only "unlinked second idiomorph" check on strong HMG hits in called
  genomes would surface them.

## Decision
None yet.

## Open
- Curate Mycotypha africana sexM (h1) as a record? (Literature: Schulz 2016.)
- Dichotomocladium elegans h5: synteny, annotation, classifier score.
- Remove probable MAT genes from the negative set; recompute the gate (approval).
- Gene-tree outgroup: use the sexP-related paralog family deliberately, or an
  HMG family chosen independently.
- Locate Bifiguratus h3 and C. minor h6.

## Files
`results/2026-10-04_unclassified_hmg/` (`NOTE.md`, `place_hmg.py`,
`hmg_outgroup_placement.tsv`, `sexP_sister_block.txt`, query proteins).
