# What marks a true pheromone-receptor call, and would ML help
Status: open (measurement only; revised after independent review, branch pr-receptor-ml-review 15a5a56)

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
### 1. Single features, headline set (20 mating, 64 other, 17 genomes = 12 species, 4 orders; `univariate_H.tsv`)
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

- CAAX proximity profile (`caax_distance_profile.tsv`, H): a strict-CAAX ORF lies within 10 kb of 12/20 mating loci and 0/64 other loci; within 50 kb, 13/20 and 7/64; within 100 kb, 14/20 and 9/64. The rule is highly specific here. It misses 8 of 20 mating loci at 10 kb and 7 of 20 at 50 kb.
- Protein shape does not separate the classes (length, helix count, PF02076 score, hydrophobicity). Intron count is weakly higher in mating loci (mean 3.3 versus 2.3), which may be an artefact of gene-model quality.
- Tree (67 complete-model loci, `tree_purity.tsv`): the nearest relative in the same order has the same label no more often than chance (mating 8/22 versus null mean 0.32). In another order, 15/22 mating loci have a mating nearest relative (null mean 0.34), but 17/45 other loci also have a mating nearest relative where the null expects far fewer (agreement for `other` 0.38, null mean 0.64). Mating receptors therefore act as hubs in the tree, close to loci of other orders, and this is not evidence of a mating-specific clade or sequence signal (200 shuffles, 22 and 45 loci).
- Flank genes, what drives the feature. The hits come from a few genes. In Tremellales the hits match the Sporidiobolales STE20 queries, in Sporidiobolales they match the Cryptococcus records (STE20a and STE20alpha), and in Agaricomycetes they match the STE20 queries. In Ustilaginales the hit lies 353 to 886 bp from every mating locus and, in U. maydis, is the PAN6 homologue (UMAG_10139, putative pantoate-beta-alanine ligase) beside rba1, not STE20 (reviewer's check of the RefSeq annotation of NC_026482.1; the STE20 hit in U. maydis is 62 kb from pra1). Cryptococcus PAN6 detecting the Ustilaginales neighbour is a cross-subphylum observation. This is the ancestral P/R-region synteny already described in the literature for Tremellales and Sporidiobolales.
- Flank genes by class: near the receptor in 8/8 Ustilaginales, 5/5 Sporidiobolales and 2/2 Tremellales mating loci, and near 2/22 Sporidiobolales and 0/26 Ustilaginales `other` loci (`n_FLANK_50kb_xo`).
- Agaricales is the fairer test. The curated flank set there (MIP1 and beta_fg, flanking the A locus) was not built from that order's P/R linkage. The flank distance AUC is 0.88 (driven by C. cinerea, where mating loci lie 87-108 kb from a flank gene and others 159 kb or more), the 50 kb count AUC is 0.56, and 1/5 mating loci has a flank hit within 50 kb.
- Position confounding (reviewer): 52 of 64 negatives lie on contigs with no flank hit at all, and "any flank hit on the same contig" alone gives AUC 0.91. Assembly fragmentation is therefore part of the signal until a position-controlled null is run.

### 2. Held-out clade comparison, headline set H, LOCO (20 mating, 64 other, 17 genomes, 12 species, 4 held-out orders)
Intervals bootstrap over genomes unless stated. Delta is the paired AUC difference to the rule.
| Model | AUC (95%) | AP | Sens at 95% spec | Delta AUC vs rule |
|---|---|---|---|---|
| Rule: strict CAAX within 10 kb | 0.80 (0.66-0.91) | 0.70 | 0.60 | 0 |
| Rule as distance to CAAX ORF | 0.79 (0.63-0.91) | 0.76 | 0.65 | -0.01 (-0.09 to 0.08) |
| Baseline, no learning: CAAX within 10 kb OR flank gene within 50 kb | 0.98 (0.96-1.00) | 0.87 | 1.00 (20/20 mating, 3/64 other flagged) | +0.18 (0.07 to 0.31) |
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
| ESM-2 8M, logistic regression | 0.38 (0.26-0.56) | 0.20 | 0.05 | -0.41 (-0.57 to -0.21) |
| ESM-2 35M, logistic regression | 0.48 (0.36-0.62) | 0.23 | 0 | -0.32 (-0.48 to -0.14) |
| ESM-2 150M, logistic regression | 0.60 (0.45-0.75) | 0.39 | 0.20 | -0.20 (-0.33 to -0.05) |
| ESM-2 650M, logistic regression | 0.53 (0.39-0.68) | 0.28 | 0.15 | -0.27 (-0.41 to -0.09) |
| ESM-2 150M, kNN | 0.65 (0.47-0.79) | 0.37 | 0.10 | -0.15 (-0.27 to -0.02) |
| ESM-2 150M plus CAAX features, first version (weighting cancelled by the scaler) | 0.71 (0.52-0.85) | 0.61 | 0.50 | -0.09 (-0.20 to 0.03) |
| ESM-2 150M plus CAAX features, corrected (blocks scaled separately, CAAX block weighted sqrt(640/6)) | 0.61 (0.40-0.77) | 0.56 | 0.50 | -0.20 (-0.35 to -0.07) |

(`cv_summary_H_LOCO.tsv`, `cv_summary_H_extra_LOCO.tsv`, `review_fixes.tsv`. Permutation p for the handcrafted and CAAX-containing models is 0.001; for sequence-only models 0.03 to 0.99, mostly not significant.)

- The review found a bug in the first ESM-2 plus CAAX model: the CAAX columns were multiplied by 3 before a StandardScaler inside the pipeline, which cancelled the weight and left 6 CAAX columns among 640 embedding columns under strong L2. In the corrected refit (`review_fixes.py`, saved embeddings, CPU, no new job) the two blocks are standardised separately and the CAAX block weighted by sqrt(640/6). The corrected AUC is 0.61, not better than the first version (0.71) and still below the rule. This test is also weak by design: logistic regression on the CAAX features alone scores 0.68 against 0.80 for the rule, so it says little about ESM-2 itself.
- Resampling unit for the lead over the rule (`review_fixes.tsv`; 17 genomes are 12 species, 5 R. toruloides strains and the near-isogenic C. neoformans pair JEC20/JEC21 count once under the species unit; the order unit has 4 orders). Delta AUC versus the rule, 95% interval:

| Model | By genome | By species | By order |
|---|---|---|---|
| Flank features only | 0.03 to 0.32 | -0.00 to 0.30 | -0.02 to 0.28 |
| All handcrafted features | 0.06 to 0.32 | 0.05 to 0.29 | 0.02 to 0.27 |
| CAAX OR flank baseline | 0.07 to 0.31 | 0.07 to 0.30 | 0.03 to 0.25 |

  The lead of the flank-only model is not robust once conspecific genomes or orders are the unit. The full handcrafted model and the simple OR baseline stay above zero under all three units. (These intervals were recomputed with a new random seed; the reviewer's values are -0.01 to 0.31 and -0.02 to 0.29 for flank only.)
- Per held-out order (`cv_per_clade_H_LOCO.tsv`; positives/negatives): rule AUC Agaricales 1.00 (5/12), Sporidiobolales 0.70 (5/22), Tremellales 0.75 (2/4), Ustilaginales 0.75 (8/26). Handcrafted logistic regression 0.93, 0.97, 1.00, 1.00. Flank features only: 0.88, 0.99, 1.00, 1.00. Without flank features: 1.00, 0.69, 0.50, 0.48. ESM-2 150M: 0.82, 0.42, 0.63, 0.58. The Tremellales fold is 6 loci from one species (two near-isogenic genomes), so its 1.00 is weak evidence.
- LOSO (secondary, leaks close homologs; `cv_summary_H_LOSO.tsv`): rule 0.80; handcrafted 0.98; ESM-2 150M 0.71 (0.59-0.80); 3-mer 0.64; first-version ESM-2 plus CAAX 0.82 (delta +0.02, -0.07 to 0.12). Sequence-only models improve when close relatives enter training.
- Sensitivity sets (LOCO): ALL (31/87, 23 genomes) rule 0.82 (0.72-0.90), handcrafted 0.92 (0.81-0.99), ESM-2 150M 0.45 (0.35-0.57). Hq (12/28) rule 0.83, handcrafted 0.96, ESM-2 150M 0.73 (0.49-0.87). Hnocc rule 0.80, handcrafted 0.98, ESM-2 150M 0.61.
- Learning curve (H, random training genomes, not split by order, `learning_curve.tsv`): handcrafted AUC on unseen genomes 0.91 with 2 training genomes, 0.96 with 4, 0.98 with 12; ESM-2 150M 0.56, 0.61, 0.68 and still rising. This curve does not show what limits ESM-2: it rises with more genomes and the data are too few to extrapolate.

## Limits
- Four held-out orders, 20 mating loci in 12 species, so intervals cannot reflect missing orders and conspecific genomes are not independent. Ustilaginales is the order with most data (8 genomes, 8 species). Pucciniomycotina rust receptors, Malasseziales, Cantharellales and most Agaricomycete orders have no labelled receptor here.
- Flank result, semantic circularity. A positive is the receptor at the curated MAT locus. The flank feature is proximity to genes that curators placed in MAT records of other lineages. The query sequences come from other orders, but the choice of which genes count as flank genes was made knowing that STE20 and PAN6 are MAT-linked in Tremellales and Sporidiobolales, two of the four test orders, so a high AUC there is close to true by definition. It confirms ancestral P/R-region synteny and is not an independent predictor of mating function. The own-species and own-order flank features give the same AUC (0.98), so no same-order query sequence is involved.
- The flank signal is one or two genes (STE20; the PAN6 homologue beside rba1 in Ustilaginales), a large part of it is contig co-occurrence, and Agaricales, the one order with an independently built flank set, gives distance AUC 0.88, count AUC 0.56 and 1/5 mating loci with a hit.
- Negatives are assumed, not tested. Grade B positives are mapped, not tested. No label was derived from a CAAX call or a PR-only pipeline call. Not checked: how each db receptor was originally chosen.
- Sequence models: the ESM-2 result is underpowered, not a demonstrated absence of signal. There are 20 positives, fewer distinct receptor sequences (5 R. toruloides alleles and two near-isogenic Cryptococcus genomes), 45 of 118 gene models lack a PF02076 hit so many embeddings come from truncated proteins, only mean pooling with linear or kNN heads was tried, and the ESM-2 plus CAAX model has a weak design. ESM-2 3B was not run.
- The widened motif feature `nT2_10kb` (AUC 0.95) was chosen with knowledge of labelled positives and is not used in models.
- Not checked: other labelled genomes outside the db, a position-controlled null for the flank feature, flank scoring of the 1,360 unverified calls, promoter or motif analysis. The conclusions concern PR calls at known mating loci.

## Recommendations/decisions
Proposals only; nothing is implemented.
1. Strict CAAX proximity is specific where it fires (12/20 mating versus 0/64 other at 10 kb) and no sequence feature tested beats it alone, but it misses 8 of 20 labelled mating loci at 10 kb.
2. On this data no sequence model (ESM-2 at four sizes, k-mers, alignment similarity, PF02076 score) reached the rule on held-out orders. With 20 positives from 4 orders the data cannot rule out a moderate sequence signal. ESM-2 should be judged again only with more labelled orders.
3. Neighbourhood is the signal that adds to the rule. The simple, unlearned baseline "CAAX within 10 kb OR flank gene within 50 kb" catches 20/20 labelled mating loci with 3/64 false positives (AUC 0.98), and the full handcrafted model adds nothing beyond it that this data can show. Read this as known ancestral MAT-region synteny (STE20, and the PAN6 homologue beside rba1), not as a discovery, and not yet as independent evidence. Any learned flank model should be compared with this baseline first.
4. Test the flank idea out of lineage before any use. (a) Build the STE20 and PAN6 query set from one lineage only (Pucciniales or Microbotryum) and score Tremellales, Ustilaginales and Sporidiobolales with it. (b) Use a null that controls for genome position: shuffle flank-gene positions within each genome, or report AUC conditional on contig length or on contigs carrying a flank hit, so that assembly fragmentation is not counted. (c) Add genomes whose receptor is independent of position (Microbotryum, Sporisorium, U. hordei, Malassezia, further Cryptococcus species).
5. Minimum labelled data before a learned classifier is worth integrating: at least 10 independent orders with 5 or more labelled mating loci each (counting species, not strains) and at least 100 tested or well-justified negatives, with receptors not chosen by CAAX. Current: 4 orders, 12 species, 20 mating loci, 64 untested negatives.
6. Safe integration, if later justified: a rank or flag column for loci already admitted by the evidence rules, rule tiers unchanged, a held-out-order check in the test suite, and no use of the score to add or remove calls.
7. Ranked next steps: (a) the out-of-lineage, position-controlled flank test above; (b) score the 1,360 unverified PR calls for flank proximity and inspect the discordant cases by hand (flank yes and CAAX no, and the reverse), reporting "CAAX OR flank" as the baseline; (c) expand the labelled set from genomes with mapped mating types and tested non-mating STE3 genes; (d) revisit ESM-2 (3B, fine-tuned or per-residue features) only after (c) reaches 10 orders.
