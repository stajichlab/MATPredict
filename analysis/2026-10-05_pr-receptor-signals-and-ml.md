# What marks a true pheromone-receptor call, and would ML help
Status: open (measurement only; no detection or database change proposed in code)

## Question
1. What labelled true-positive and true-negative receptor loci exist, and how good are the labels?
2. Which measurable properties separate mating receptors from non-mating STE3-like copies? Is proximity to a CAAX ORF the dominant signal?
3. Would a classifier on protein sequences, including ESM-2 embeddings, beat the current rule (STE3-like locus with a strict CAAX ORF within 10 kb) on held-out clades?

## Data and label provenance
Files: `results/2026-10-05_pr_receptor_ml/labelled_loci.tsv` (one row per STE3-like locus, with provenance, grade and features), `labelled.faa`, `protein_features.tsv`.

- Loci: miniprot of the 1,094 STE3 queries of PR #33 (`scan_genome.py`) against 23 genomes, merged as in PR #33. 118 loci: 31 mating, 87 other.
- Mating label, no CAAX evidence used. Agaricomycete panel genomes: overlap with the curated B-locus receptor cluster (PR #33 `panel_loci.tsv`). All other genomes: the locus whose best curated-reference hit is the genome's own db record receptor at 95 percent identity or more. Every other STE3-like locus is `other`.
- Grades of the mating label: A genetically tested in that species (U. maydis pra1, PMID 1310895; C. neoformans STE3alpha and STE3a, PMID 12455690; S. commune bar3/bbr2, PMID 7489716): 5 loci. B mapped to the mating locus, not tested in that species (7 Ustilaginales a-locus records, 5 Rhodotorula P/R records, C. cinerea B43): 15 loci. C inferred or putative (Trametes, Grifola, Russula, Heterobasidion records chosen from CAAX positional evidence; two Wallemia putative records): 11 loci. D `other`: assumed non-mating, never tested.
- Circularity guard. Grade C genomes are not used in the headline set because their positives were chosen with CAAX evidence, and the PR-only calls of the pipeline are never used as labels. Headline set H = grades A and B: 20 mating, 64 other, 17 genomes, 4 orders (Agaricales, Sporidiobolales, Tremellales, Ustilaginales). Sensitivity sets: ALL (31/87, 23 genomes, 7 orders), Hq (complete gene models, 12/28, 13 genomes), Hnocc (H without the C. cinerea `other` copies, 20/59).
- Negatives are weak. C. cinerea carries about 14 receptor genes, so some of its 5 `other` copies may be mating. No negative has been tested.

## Method
- Features per locus (all measured, none assumed): strict-CAAX count and nearest distance (`d_T`, 10/20/50 kb), precursor homology (tblastn of curated precursors), nearest HD-class gene and nearest conserved flank gene (STE20, RPO41, UAP1, FAO1, NOG2, PAN6, RPL39, MIP1, beta_fg; miniprot of the db record proteins), protein length, intron count, transmembrane helices (Kyte-Doolittle), PF02076 score and coverage (pyhmmer, Pfam PF02076.23).
- Leakage controls. The curated queries of the genome's own order are removed before the precursor, HD and flank features are computed (`recompute_xo.py`), so a new order is treated as having no records. Gene models for protein features come from queries without the genome's own record (`reprotein.py`). `best_ident`, `best_ref_ident`, `n_queries` and `aln_score` in the TSV are leaky and are not used.
- Evaluation: leave-one-order-out (LOCO, headline) and leave-one-genome-out (LOSO, secondary, leaks close homologs inside an order). Out-of-fold scores are pooled. Intervals: bootstrap over genomes (1,000 or 2,000 draws); permutation of labels within genome for the null; paired bootstrap of the AUC difference to the rule. Sensitivity at 95 percent specificity is the best true-positive rate with false-positive rate at or below 0.05.
- Models: rule (strict CAAX within 10 kb); logistic regression, gradient boosting and random forest on handcrafted features; ESM-2 8M, 35M, 150M and 650M mean-pooled embeddings with logistic regression or kNN; 2-mer and 3-mer logistic regression; nearest-neighbour Smith-Waterman similarity (`alignNN`). Embeddings and models ran on CPU (HPCC `short`).
- Tree: FastTree of the 67 complete-model loci plus S. cerevisiae Ste3; nearest-relative label agreement.

## Results
### 1. Single features, headline set (20 mating, 64 other, 17 genomes; `univariate_H.tsv`)
| Feature | AUC (95% bootstrap over genomes) | Within-genome AUC |
|---|---|---|
| Strict CAAX within 10 kb (the rule) | 0.80 (0.67-0.91) | 0.77 |
| Distance to nearest strict-CAAX ORF (closer = mating) | 0.79 (0.63-0.91) | 0.78 |
| Precursor homology within 10 kb, own order removed | 0.56 (0.49-0.66) | 0.58 |
| Conserved flank genes within 50 kb, own order removed | 0.87 (0.75-0.99) | 0.94 |
| Distance to nearest flank gene (closer = mating) | 0.98 (0.96-1.00) | 0.99 |
| HD-class gene within 20 kb | 0.52 (0.48-0.58) | 0.53 |
| Protein length | 0.50 (0.37-0.62) | 0.47 |
| Intron count | 0.70 (0.56-0.84) | 0.71 |
| Transmembrane helices | 0.56 (0.41-0.71) | 0.49 |
| PF02076 score | 0.50 (0.38-0.63) | 0.48 |
| Max similarity to a locus of another genome, same order | 0.38 (0.24-0.56) | 0.34 |

- CAAX proximity profile (`caax_distance_profile.tsv`, H): a strict-CAAX ORF lies within 10 kb of 12/20 mating loci and 0/64 other loci; within 50 kb, 13/20 and 7/64; within 100 kb, 14/20 and 9/64. The rule is highly specific here and misses 8 of 20 mating loci at any distance up to 50 kb.
- Protein shape does not separate the classes (length, helix count, PF02076 score, hydrophobicity). Intron count is weakly higher in mating loci (mean 3.3 versus 2.3), which may be an artefact of gene-model quality.
- Tree (67 loci, `tree_purity.tsv`): the nearest relative in the same order has the same label no more often than chance (mating 8/22 versus null mean 0.32). The nearest relative in another order is mating for 15/22 mating loci (null mean 0.34, 95th percentile 0.64), a weak hint that mating receptors of different orders sit closer to each other than expected; the interval is not tight (200 shuffles, 22 loci).
- Flank genes (STE20, RPO41, MIP1 and others) are near the mating receptor in all 8 Ustilaginales, 5 Sporidiobolales and 2 Tremellales mating loci, and near 2 of 22 Sporidiobolales and none of 26 Ustilaginales `other` loci (`labelled_loci.tsv`, `n_FLANK_50kb_xo`). In U. maydis STE20 lies 62 kb from pra1 (checked directly with miniprot). This is conserved synteny of the P/R region, not proximity to a CAAX ORF.

### 2. Held-out clade comparison, headline set H, LOCO (20 mating, 64 other, 17 genomes, 4 held-out orders)
| Model | AUC (95%) | AP | Sens at 95% spec | Delta AUC vs rule (95%) |
|---|---|---|---|---|
| Rule: strict CAAX within 10 kb | 0.80 (0.66-0.91) | 0.70 | 0.60 | 0 |
| Rule as distance to CAAX ORF | 0.79 (0.63-0.91) | 0.76 | 0.65 | -0.01 (-0.09 to 0.08) |
| Logistic regression, CAAX and precursor features | 0.68 (0.51-0.83) | 0.61 | 0.50 | -0.11 (-0.19 to -0.03) |
| Logistic regression, all handcrafted features | 0.97 (0.95-1.00) | 0.86 | 0.85 | +0.18 (0.06 to 0.32) |
| Logistic regression, no CAAX or precursor features | 0.95 (0.89-1.00) | 0.80 | 0.80 | +0.16 (0.00 to 0.32) |
| Logistic regression, flank features only (2 features) | 0.96 (0.92-1.00) | 0.81 | 0.80 | +0.16 (0.01 to 0.32) |
| Logistic regression, everything except flank features | 0.63 (0.47-0.79) | 0.56 | 0.45 | -0.17 (-0.28 to -0.07) |
| Logistic regression, HD features only | 0.44 (0.34-0.52) | 0.22 | 0 | -0.36 (-0.48 to -0.23) |
| Random forest, all handcrafted | 0.97 (0.94-1.00) | 0.87 | 1.00 | +0.18 (0.06 to 0.31) |
| Gradient boosting, all handcrafted | 0.96 (0.91-1.00) | 0.84 | 0.75 | +0.16 (0.02 to 0.32) |
| 3-mer composition, logistic regression | 0.53 (0.41-0.68) | 0.26 | 0 | -0.26 (-0.41 to -0.10) |
| Smith-Waterman nearest neighbour | 0.61 (0.50-0.75) | 0.30 | 0 | -0.19 (-0.38 to 0.03) |
| PF02076 score alone | 0.50 (0.37-0.63) | 0.22 | 0 | -0.30 (-0.48 to -0.10) |
| ESM-2 35M, logistic regression | 0.48 (0.36-0.62) | 0.23 | 0 | -0.32 (-0.47 to -0.12) |
| ESM-2 150M, logistic regression | 0.60 (0.45-0.75) | 0.39 | 0.20 | -0.20 (-0.33 to -0.05) |
| ESM-2 650M, logistic regression | 0.53 (0.39-0.68) | 0.28 | 0.15 | -0.27 (-0.41 to -0.09) |
| ESM-2 150M, kNN | 0.65 (0.47-0.79) | 0.37 | 0.10 | -0.15 (-0.27 to -0.02) |
| ESM-2 150M plus CAAX features | 0.71 (0.52-0.85) | 0.61 | 0.50 | -0.09 (-0.20 to 0.03) |

(`cv_summary_H_LOCO.tsv`, `cv_summary_H_extra_LOCO.tsv`; ESM-2 8M also in the table files, AUC 0.38. Permutation p for the models above 0.9: 0.001; for the sequence-only models p is 0.03 to 0.9, mostly not significant.)

- Per held-out order (`cv_per_clade_H_LOCO.tsv`; positives/negatives): rule AUC Agaricales 1.00 (5/12), Sporidiobolales 0.70 (5/22), Tremellales 0.75 (2/4), Ustilaginales 0.75 (8/26). Handcrafted logistic regression 0.93, 0.97, 1.00, 1.00. Flank features only: 0.88, 0.99, 1.00, 1.00. Without flank features: 1.00, 0.69, 0.50, 0.48. ESM-2 150M: 0.82, 0.42, 0.63, 0.58.
- LOSO (secondary, leaks close homologs; `cv_summary_H_LOSO.tsv`): rule 0.80; handcrafted 0.98; ESM-2 150M 0.71 (0.59-0.80); 3-mer 0.64; ESM-2 plus CAAX 0.82 (0.72-0.91, delta to rule +0.02, -0.07 to 0.12). Sequence-only models improve under LOSO as close relatives enter training, and still stay below the rule.
- Sensitivity sets (LOCO): ALL (31/87, 23 genomes) rule 0.82 (0.72-0.90), handcrafted 0.92 (0.81-0.99), ESM-2 150M 0.45 (0.35-0.57). Hq (12/28) rule 0.83, handcrafted 0.96, ESM-2 150M 0.73 (0.49-0.87). Hnocc rule 0.80, handcrafted 0.98, ESM-2 150M 0.61.
- Learning curve (H, random training genomes, `learning_curve.tsv`): handcrafted AUC on unseen genomes 0.91 with 2 training genomes, 0.96 with 4, 0.98 with 12; ESM-2 150M 0.56, 0.61, 0.68; the rule has no training and stays at about 0.79. Training size is not what limits ESM-2.

## Limits
- Four held-out orders, 20 mating loci, 7 independent genomes in the order with most data. Intervals bootstrap over genomes but cannot reflect orders that are missing. Pucciniomycotina rust receptors, Malasseziales, Cantharellales and most Agaricomycete orders have no labelled receptor here.
- The handcrafted gain rests on two features, proximity to conserved flank genes. The flank queries come from curated records of the same kind of locus (Tremellales STE20 and RPO41 genes, Sporidiobolales STE20, Agaricomycete MIP1 and beta_fg). The own order is removed, so no query comes from the same order, but ancestral synteny between orders (Tremellales, Sporidiobolales, Ustilaginales) makes the transfer easy and it is weaker in Agaricales (flank-only AUC 0.88 there; 1 of 5 mating and 1 of 12 other loci have a flank hit within 50 kb) and cannot be tested in orders without such a record. This is a measured synteny signal, not an independent validation.
- Negatives are assumed, not tested. Grade B positives are mapped, not tested. No label was derived from a CAAX call or from a PR-only pipeline call, but the rule is still scored against labels from genomes whose records were curated by authors who may have used pheromone-linked evidence (not checked per record).
- Gene models are miniprot models and many are truncated (45 of 118 have no PF02076 hit); protein features and embeddings are therefore noisy. Hq restricts to complete models and shows the same ranking.
- The motif-widened `nT2_10kb` (AUC 0.95) was chosen with knowledge of labelled positives and is not used in models.
- ESM-2 3B was not run (the GPU job was cancelled before it started); 8M, 35M, 150M and 650M were run on CPU. Only mean pooling and linear or kNN heads were tried. No fine-tuning (too few labels).
- Not checked: other labelled genomes outside the db (literature), per-record evidence of how each db receptor was chosen, other classifiers on the tree, a promoter or motif analysis.
- The conclusions concern PR calls at known mating loci. They do not test the 1,360 unverified calls directly.

## Recommendations/decisions
Proposals only; nothing is implemented.
1. Proximity to a strict-CAAX ORF is a strong specific signal where it fires (12/20 versus 0/64 at 10 kb) and no tested sequence feature beats it alone. It is not the whole story: 40 percent of labelled mating loci have none within 10 kb.
2. Sequence-based ML does not improve detection on this data. ESM-2 (four sizes), k-mers and alignment similarity all fall below the rule on held-out orders; the best ESM-2 interval excludes parity in the wrong direction. Adding ESM-2 to CAAX features does not beat the rule (delta -0.09, -0.20 to 0.03). Receptor protein sequence carries little information about mating function here, consistent with the unresolved PF02076 tree.
3. The most promising signal is neighbourhood, not sequence: conserved MAT-region flank genes (STE20, RPO41 and relatives) next to the receptor. It beats the rule on held-out orders (delta AUC +0.16, 0.01 to 0.32) but with 4 orders and one source of queries it is a lead to test, not a result to deploy.
4. Minimum data before a learned classifier is worth integrating: at least 10 independent orders with 5 or more labelled mating loci each and at least 100 tested or well-justified negatives, with receptors not chosen by CAAX. Current: 4 orders, 20 mating, 64 untested negatives.
5. Safe integration, if later justified: output a score and rank for loci already admitted by the evidence rules, with the rule tiers unchanged, a held-out-order check in the test suite, and no use of the score to add or remove calls.
6. Ranked next steps: (a) add labelled loci from genomes with mapped mating types (Microbotryum, Sporisorium, Ustilago hordei, Malassezia, Cryptococcus species, Pucciniales STE3 loci, Agaricomycete B-locus genomes beyond the six used) and tested non-mating STE3 genes; (b) curate the flank-gene set for each lineage independently of the receptors and measure a synteny score on the 1,360 unverified calls in orders where the flank is known; (c) report the flank score beside the CAAX flag in `loci.tsv` as an extra confirmation column; (d) revisit ESM-2 only after (a) gives 10 or more orders.
