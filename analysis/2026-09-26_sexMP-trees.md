# sexM/sexP candidate phylogenies
Status: decided (tree is supporting evidence, not the classifier)

## Question
Do candidate HMG genes group with curated sexM or sexP, and can a tree support
tier-2 curation decisions?

## Data and code version
Alignment: hmmalign to Pfam PF00505 (69 columns). Candidates from detection
reports (early-diverging scan, then the classifier scans).

## Method
1. IQ-TREE on 878 tips (`results/2026-09-26_sexMP_phylogeny/`).
2. FastTree redraw with classifier labels, 868 tips (`results/2026-09-27_sexMP_fasttree/`).
3. Trimmed final ML tree, 181 tips, IQ-TREE (1000 UFBoot) and RAxML-NG (500
   bootstraps), rooted on the F. graminearum MAT1-2-1 (`results/2026-09-27_sexMP_final_ml/`).

## Results
- 878-tip IQ-TREE: sexP clade UFBoot 99; sexM refs not monophyletic (7/9, UFBoot 57).
- FastTree: label/clade agreement 148/195 before, 204/206 after the classifier.
- Final trimmed tree: sexP 62 UFBoot / 6 bootstrap; sexM (all 9 refs, first
  time monophyletic) 74 / 0. Label/clade agreement 94/97. No
  Lichtheimiaceae/Syncephalastraceae copy reaches UFBoot >= 95.

## Limits
One 69-column domain cannot give strong support even for the references.

## Curator decisions
Made: the tree is supporting evidence; the UFBoot>=95 curation rule was
replaced by classifier margin >25 bits plus gene order (2026-09-27).

## Files
`results/2026-09-26_sexMP_phylogeny/tree_hmm.pdf`,
`results/2026-09-27_sexMP_fasttree/tree_ft.pdf`,
`results/2026-09-27_sexMP_final_ml/iq.pdf`, each folder's `tip_names.tsv` and `NOTE.md`.
