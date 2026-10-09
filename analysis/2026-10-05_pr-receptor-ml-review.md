# Review: PR receptor signals and ML (branch pr-receptor-ml-investigation, fb2ab09)
Status: independent review of `analysis/2026-10-05_pr-receptor-signals-and-ml.md`. All numbers below were recomputed from the committed tables in `results/2026-10-05_pr_receptor_ml/`. No HPCC jobs were run and nothing in the study was changed.

## Verdict
- The counts, AUCs and intervals in the report match the committed tables.
- The main conclusion holds for this dataset: on held-out orders, protein sequence models do not beat the CAAX rule, and the conserved-neighbour (flank) feature does.
- The flank result needs softer wording in three places: what the feature actually is, how independent it is from the labels, and how firm its lead over the rule is.
- The ESM-2 "no" is a "no evidence on 20 positives", not evidence that sequence carries no signal.

## 1. Label leakage and circularity
- No direct leakage. Labels come from receptor identity to the genome's own record (at least 95 percent) or from overlap with the panel B-locus cluster. Flank genes are never used to make a label (`make_labels.py`).
- The own-order exclusion works (`recompute_xo.py`). Every taxid that supplies a FLANK query is in `taxid_order.tsv`, and so is every taxid in the precursor set. The own-species and own-order flank features give the same AUC (0.983, within-genome 0.992), so all of the signal comes from other orders.
- There is a weaker, semantic circularity. A positive means "the receptor at the curated MAT locus". The flank feature means "next to genes that curators put in MAT records of other lineages". The query sequences come from other orders, but the choice of which genes count as flank genes was made knowing that STE20 and PAN6 are MAT-linked in Tremellales and Sporidiobolales, which are two of the four test orders. Removing the query sequences does not remove that prior knowledge. In those lineages a high AUC is close to true by definition: it confirms ancestral P/R-region synteny that is already in the literature, and it is not an independent predictor of mating function.
- The feature is really one or two genes, judging from the hit counts per locus:
  - Tremellales: 5 hits, which matches the 5 Sporidiobolales STE20 records.
  - Sporidiobolales: 2 hits, which matches the 2 C. neoformans records (STE20a and STE20alpha are the likely genes).
  - Agaricomycetes: 7 hits, which matches all 7 STE20 queries.
  - Ustilaginales: 2 hits, at 353 to 886 bp from every mating locus. For U. maydis I checked the RefSeq annotation of NC_026482.1. The hit 803 bp from the locus is UMAG_10139, a putative pantoate-beta-alanine ligase (the PAN6 homologue) next to rba1, not STE20. The report's description ("STE20, RPO41, MIP1 and others"; STE20 at 62 kb) does not name the gene that drives the Ustilaginales signal.
- A fairer estimate of the flank AUC is in Agaricales, the only order where the curated flank set (MIP1 and beta_fg, which flank the A locus) was not built from that order's own P/R linkage:
  - distance AUC 0.88, which comes from ranking at about 100 kb in C. cinerea (mating 87 to 108 kb, others 159 kb or more);
  - AUC 0.56 when counting flank hits within 50 kb;
  - 1 of 5 mating loci has a flank hit within 50 kb.
  - Grade C genomes are not usable for this check, because their labels came from CAAX evidence.

## 2. Genome confounding and LOCO with 4 orders
- Genome-level differences do not explain the flank AUC. All 17 genomes in set H carry both classes, and the within-genome AUC is 0.99.
- Much of the AUC comes from contig co-occurrence. 52 of 64 negatives sit on contigs with no flank hit at all, and "any flank hit on the same contig" alone gives AUC 0.91.
- The genomes are not independent. The 17 genomes are only 12 species: 5 strains of R. toruloides, and the near-isogenic pair C. neoformans JEC20 and JEC21. The Tremellales fold is 6 loci from one species, so its per-fold AUC of 1.00 tells us little.
- The lead of the 2-feature flank model over the rule (dAUC +0.16) depends on the bootstrap unit:

| Bootstrap unit | dAUC vs rule (95% interval) |
|---|---|
| genome (as reported) | 0.02 to 0.32 |
| species | -0.01 to 0.31 |
| order | -0.02 to 0.29 |

  - So "beats the rule" is not robust once conspecific genomes are not counted as independent.
  - The full handcrafted model (`LR_all`) stays above zero under all three units: 0.05 to 0.32, 0.05 to 0.29 and 0.02 to 0.27.

## 3. Was all same-order information removed?
- The query sequences were removed correctly.
- Same-order information is still present in two other ways:
  - the list of flank genes itself (see section 1);
  - the cross-testing between Tremellales and Sporidiobolales, where each order's STE20 is the query for the other. These orders are in different subphyla, so this is genuine cross-order transfer, but it is the case most likely to succeed.
- In Ustilaginales the signal comes from Cryptococcus PAN6, which is a real cross-subphylum observation and worth reporting.

## 4. ESM-2 fairness
- The setup is reasonable:
  - mean pooling of the last layer, with special tokens excluded;
  - standardised features, L2 logistic regression with C = 0.01 and balanced class weights;
  - a cosine kNN head as an alternative.
- Mixing score scales across folds does not explain the low ESM-2 numbers. The ESM-2 150M AUC is 0.60 pooled and 0.60 when scores are ranked within each order first.
- The ESM-2 "no" is underpowered:
  - There are 20 positives, but fewer distinct receptor sequences, because 5 R. toruloides alleles and 2 near-isogenic Cryptococcus genomes are among them.
  - 45 of 118 gene models lack a PF02076 hit, so many embeddings come from truncated proteins.
  - The learning curve still rises (0.56 to 0.68 from 2 to 12 training genomes), and the LOSO AUC is 0.71. The sentence "Training size is not what limits ESM-2" is not supported. That learning curve also splits by random genome, not by order.
- The "ESM-2 plus CAAX" model is handicapped by its design. The `*3.0` weighting of the CAAX columns has no effect, because `StandardScaler` runs after `hstack` inside the pipeline. The 6 CAAX features are therefore drowned among 640 embedding dimensions under strong L2. Logistic regression on the CAAX features alone already scores only 0.68 against 0.80 for the rule, so "adding ESM-2 does not beat the rule" says little about ESM-2.
- The tree "hint" is ambiguous. The nearest other-order relative of an `other` locus is also unusually often a mating locus (17 of 45 have the same label, against a null mean of 0.64). Mating receptors look like hubs in the tree rather than a clade, so the 15 of 22 figure is not evidence of a mating-specific sequence signal.

## 5. Recomputation and mismatches
- Confirmed:
  - 118 loci (31 mating, 87 other) in 23 genomes;
  - set H: 84 loci (20 mating, 64 other) in 17 genomes and 4 orders;
  - grades: 5 A, 15 B, 11 C;
  - rule AUC 0.80 (0.67 to 0.91);
  - flank distance AUC 0.983 (0.956 to 1.00), within-genome 0.99;
  - flank-only logistic regression under LOCO: 0.957 when I refit it;
  - the LOCO, LOSO, ALL, Hq and Hnocc table values quoted in the report;
  - the per-order AUCs;
  - the flank-hit counts (2 of 22 Sporidiobolales and 0 of 26 Ustilaginales `other` loci).
- Wrong:
  - "misses 8 of 20 mating loci at any distance up to 50 kb". Within 50 kb, 13 of 20 have a CAAX ORF, so the rule misses 7 at 50 kb and 8 at 10 kb.
  - "7 independent genomes in the order with most data". Ustilaginales has 8 genomes, from 8 species.
- Not shown in the report: "CAAX within 10 kb OR flank gene within 50 kb", a simple rule with no learning, catches 20 of 20 mating loci with 3 of 64 false positives (AUC 0.98). That is the fair baseline for any learned flank model.

## 6. Wording to soften
- "The most promising signal is neighbourhood ... beats the rule (+0.16, 0.01 to 0.32)": add that the interval crosses zero under a species or order bootstrap.
- Name the genes that drive the signal: STE20, and in Ustilaginales a PAN6 homologue. Present the result as known ancestral MAT synteny, not a discovery.
- "Receptor protein sequence carries little information about mating function": change to "no sequence model reached the rule with 20 positives from 4 orders; the data cannot rule out a moderate signal".
- "Training size is not what limits ESM-2": remove it.
- "the best ESM-2 interval excludes parity in the wrong direction": true for ESM-2 150M LR alone, not for ESM-2 plus CAAX, whose interval crosses zero.

## What would settle the flank-synteny question
1. Define the flank set without the test lineages. Build the STE20 and PAN6 query set from one lineage only (for example Pucciniales or Microbotryum), and score Tremellales, Ustilaginales and Sporidiobolales with it.
2. Test genomes where the label is independent of position: more grade-A species with a genetically tested receptor, including Microbotryum, Sporisorium, U. hordei and Malassezia.
3. Use a null that controls for genome position. Shuffle flank-gene positions within each genome, or report the AUC conditional on contig length or on being on a contig with a flank hit, so that assembly fragmentation is not counted as signal.
4. Score the 1,360 unverified PR calls for flank proximity and look at the discordant cases (flank yes and CAAX no, and the reverse) by hand.
5. Report "CAAX OR flank" as the baseline before any learned model.
