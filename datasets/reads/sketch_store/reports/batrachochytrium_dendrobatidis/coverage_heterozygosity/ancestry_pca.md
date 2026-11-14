# Bd isolate ordination (PCA) from SNP allele fractions (2026-11-14)

Status: exploratory. Script `scripts/bd_snp_pca.py` (test in `tests/test_bd_coverage.py`); tables `tables/pca_scores.tsv`, `tables/pca_variance.tsv`.
Input: the 69 per-isolate SNP tables against CLFT044. Method: biallelic SNPs on the 14 scaffolds of 500 kb or more, DP >= 10 and within 0.5x to 2x of the
isolate's median depth, every 10th position, sites carried by 3 or more isolates and by at most 90% of isolates (15,804 sites). Entry = alternative-allele
fraction; an isolate with no SNP row at a site gets 0. Mean-centred, SVD. This is allele-fraction space, not a tree.

## Result
- PC1 holds 47.3% of the variance (PC2 to PC4: 4.3, 2.8, 2.4%). It is the CLFT-type against JEL-type axis: CLFT-type isolates at +26 to +52, JEL-type mostly at -23 to -9.
- PC4 separates two isolates from all others: KRBOOR_331 (South Korea) and the Swiss isolate, both at about -31 (all others between -4 and +4). Both have high depth (171x and 81x), so this is not a depth effect. This is the third-lineage-like signal behind their low JEL423 recovery in the parent-of-origin analysis.
- CLFT160 (HET-rich; 26x) lies at PC1 +16.5, between the class centres (JEL-type median -17.1, CLFT-type median +51.4), about halfway, and at PC4 +2.2, so not on the Swiss/Korean axis. JEL-type isolates at similar depth (22 to 30x) lie at -21.8 to -7.1, so the position is not a depth effect. This agrees with a mixture of the CLFT-type and JEL-type lineages and does not fit the third lineage.
- CLFT039 (HET-rich; 15x) lies at PC1 +21.1 and PC4 -1.05. This position cannot be read: see the depth effect below.

## The depth effect (a limit of this ordination)
A missing SNP row is scored as the reference allele, so low depth pulls an isolate toward the middle of PC1. Within the JEL-type class, depth and PC1 correlate with Spearman rho = -0.78 (n = 52). JEL-type isolates of 17x or less have a median PC1 of -1.2 (range -8.6 to +23.6) against -17.6 at 25 to 35x and -20.8 at 60x or more. CLFT039 (15x) falls inside the range of those low-depth isolates, so the ordination does not separate it from them; six other isolates of 13 to 16x sit between the class centres for the same reason and are not evidence of mixture.

## What this does not show
- A tree, a formal ancestry model (no admixture proportions), or a third-lineage reference. The statement "the lineage mixture fits CLFT160 and not the Swiss/Korean third-lineage signal" rests on one isolate at adequate depth.
- Anything about CLFT039 beyond the parent-of-origin windows (see `report.md`), which are not depth-corrected either.
