# sexM vs sexP profile HMMs: a discrimination test

2026-09-26. Script `hmm_loo.py`; output `results.txt`, `loo_margins.tsv`,
`disputed_scores.tsv`. HMMER 3.4, MAFFT 7.505 (cluster installs; MAFFT needs
`MAFFT_BINARIES=/opt/linux/rocky/8.x/x86_64/pkgs/mafft/7.505/libexec/mafft`).

## Design, and how it avoids the 2026-09-21 failure

The Saccharomyces a1/alpha1 HMM failed on (1) near-identical training
sequences and (2) searching 6-frame ORFs. Here:

- Every HMM scores **modelled proteins** (miniprot/tblastn models from the
  phylogeny run, funannotate proteins for Zygo). So this tests discrimination
  only, not sensitivity on genomic DNA.
- **Truth set, 38 proteins:** 15 curated references (9 sexM, 6 sexP) plus the
  single HMG-box protein annotated inside each of the 23 Zygo truth loci
  (7 Minus = sexM, 16 Plus = sexP). Each Zygo locus has exactly one HMG-box
  protein (PF00505, E <= 1e-3).
- **Diversity** (pairwise identity over co-aligned columns, min / median / max):
  sexM 16 seqs, 16 unique, 17.1 / 29.5 / 100; sexP 22 seqs, 18 unique,
  24.9 / 31.8 / 100. Diverse, unlike the Saccharomyces set (78-99%).
  Caveat: 12 genera only, all Mucoromycota; Cunninghamella gives 9 of 38.
- **Leave-one-genus-out.** sexP was also trained with the 121 sexP-clade tree
  members (UFBoot 99) minus the held-out genus and minus every disputed protein
  (`_aug` arms). sexM was never augmented (its clade has UFBoot 57).
- Baseline: blastp of each protein against the same training proteins, best
  bitscore per type. This approximates detect's best-bitscore rule on
  proteins; it is not detect's tblastn on genomes.

## Results

| arm | correct | min margin on correct side | median |
|---|---:|---:|---:|
| full-length HMM | 38/38 | 36.6 | 121.0 |
| HMG-box HMM | 38/38 | 14.8 | 60.5 |
| full-length HMM, sexP augmented | 38/38 | 33.4 | 124.9 |
| HMG-box HMM, sexP augmented | 38/38 | 16.2 | 57.5 |
| blastp baseline | 38/38 | **5.1** | 124.8 |

All methods separate sexM from sexP when given full modelled proteins. The
HMMs give a wider worst-case margin than blastp (36.6 vs 5.1 bits), the same
pattern as 2026-09-21: better discrimination when the model has a protein.

## Disputed sets (models on all 38; augmented models exclude disputed proteins)

| set | n | full HMM | HMG HMM | blastp |
|---|---:|---|---|---|
| Minus calls whose HMG box is in the sexP clade | 13 | 13 sexP (margin 13.0-191.8) | 13 sexP | 13 sexP |
| Rhizopus sexM-only calls labelled Plus | 33 | 33 sexM (-116.8 to -150.4) | 33 sexM | 33 sexM |
| Lichtheimiaceae / Syncephalastraceae sexP-clade loci | 3 | 2 sexP, 1 sexM | 2 sexP, 1 sexM | 3 sexP |

- The 13 "Minus" calls are sexP by every method: they are Plus loci labelled
  Minus. The 33 Rhizopus "Plus" calls are sexM by every method: Minus loci
  labelled Plus (the btbA vote). Both are labelling errors in detect, not
  hard discrimination cases.
- Syncephalastrum racemosum (called) and Rhizomucor pusillus (withheld): sexP
  by every method, margins 23-125. Lichtheimia ornata (withheld): unresolved
  (full HMM -8.9, HMG HMM -5.7, blastp +0.8).
- Tree vs HMM on 175 tree-placed candidates: 175/175 agree (sexM clade 51,
  sexP clade 124).

## How definitive is tree placement?

For sexP, membership in the UFBoot-99 clade is strong evidence, and every
independent method agrees with it here. For sexM it is not definitive: the
references are not monophyletic (best clade 7 of 9, UFBoot 57), so "not in the
sexP clade" does not mean sexM, and individual placements inside either clade
rest on a 69-column domain with low internal support. Use the tree as
corroboration, not as the assigner.

## Is an HMM step worth adding to detect?

Not as a fix for the mislabels: blastp on the same modelled proteins already
gets all 49 disputed-set proteins the same way. The errors come from (1) btbA
voting and (2) only one of sexM/sexP being modelled, so the decision rests on a
thin first-pass tblastn margin. Modelling both genes and comparing the
proteins fixes both. An HMM comparison of the two models is a reasonable
tie-breaker, because its worst-case margin is 7x wider than blastp's, but it
would only matter for close cases such as L. ornata. Untested: genomes outside
Mucoromycota, and loci where the model is a fragment.
