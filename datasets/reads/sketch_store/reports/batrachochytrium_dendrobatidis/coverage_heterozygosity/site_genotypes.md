# Bd depth-controlled genotypes at the divergent sites (2026-11-14)

Status: exploratory. Scripts `scripts/bd_site_genotype.slurm`, `bd_site_allele_counts.py`, `bd_site_summary.py` (tests in `tests/test_bd_coverage.py`); job array 32093759.
Tables: `tables/site_genotype_summary.tsv`, `tables/site_genotype_by_sequence.tsv`.

## Why
The SNP-table analyses (parent of origin, ordination) score a missing SNP row as the reference allele, so low depth biases them (`ancestry_pca.md`: depth and PC1
correlate at rho -0.78 in the JEL-type class). Here every isolate is re-mapped and downsampled to the same depth, and reads are counted at the same sites.

## Method
Nine isolates were chosen: the two HET-rich isolates (CLFT160, CLFT039), CLFT165 (CLFT-type with an intermediate depth profile), CLFT044 (CLFT-type control, the
reference lineage) and five JEL-type isolates (CLFT130, CLFT120, CLFT088, CLFT166, CLFT156; 13 to 29x). For each: read-count check against ENA (all matched), bwa-mem2 to CLFT044,
duplicate removal, random downsampling to 12x of the median depth, `samtools mpileup` (MAPQ and base quality 20, no BAQ) at the 172,706 single-base differences between
JEL423 and CLFT044. Reads with the CLFT044 base and reads with the JEL423 base are counted; the JEL fraction f = n_jel / (n_clft + n_jel) is taken at sites with at least 5 reads.
Classes: clft f <= 0.2, het 0.2 < f < 0.8, jel f >= 0.8.

## Result (12x for every isolate)
| isolate | class | sites with >= 5 reads | CLFT-only | heterozygous | JEL-only |
|---|---|---|---|---|---|
| CLFT044 | CLFT-type (reference lineage) | 136,328 | 0.890 | 0.105 | 0.005 |
| CLFT165 | CLFT-type | 119,222 | 0.670 | 0.244 | 0.086 |
| CLFT160 | HET-rich | 122,521 | 0.414 | 0.468 | 0.118 |
| CLFT039 | HET-rich | 122,355 | 0.408 | 0.440 | 0.152 |
| CLFT130, CLFT120, CLFT088, CLFT166, CLFT156 | JEL-type | 94,652 to 104,177 | 0.186 to 0.209 | 0.150 to 0.189 | 0.602 to 0.662 |
- The two HET-rich isolates are separated from both control groups at equal depth. Their JEL-only share (0.12 and 0.15) is a fifth of the JEL-type value (0.60 to 0.66), their heterozygous share (0.44 and 0.47) is 2.3 to 3.1 times higher than the JEL-type controls (0.15 to 0.19) and four times that of CLFT044 (0.105).
- CLFT039 (15x; the isolate the ordination could not place): per-scaffold heterozygous share on the 14 scaffolds of 500 kb or more is 0.36 to 0.56 (median 0.46) on every one of them. No control isolate exceeds 0.30 on any scaffold (JEL-type maximum 0.29, CLFT044 0.30). So CLFT039 is heterozygous for both lineages across the whole genome and does not overlap the controls. The depth confound of the ordination is removed here: its position is not a low-depth effect.
- CLFT160: median 0.49, but patchy (0.08 to 0.73 across the 14 scaffolds); low-heterozygosity scaffolds fit the two CLFT-like scaffolds seen with the parent-of-origin windows (`aneuploidy.md`).
- CLFT165: median 0.24 (0.07 to 0.53), mostly CLFT-only with an excess of heterozygous and JEL-only sites compared with CLFT044 (JEL-only 0.086 against 0.005). That is a partial mixture; I did not examine which regions.

## Limits
- JEL-type isolates show 19 to 21% CLFT-only sites and 60 to 66% JEL-only sites, not 100% JEL: shared alleles between the two lineages, mapping to a CLFT044 reference and the site choice (sites from one alignment) all contribute, so the shares are relative, not mixture proportions.
- Only nine isolates; two HET-rich isolates and one intermediate isolate. No isolate of known hybrid status, no third-lineage reference, no admixture model. A third lineage that shares alleles with both references would give heterozygous sites too; the ordination put the HET-rich isolates off the Swiss/Korean axis (`ancestry_pca.md`), but at 15x for CLFT039.
- Downsampling to 12x gives binomial noise at a site (the heterozygous class is wide, 0.2 to 0.8): a site with true f of 0.5 lands in the het class with high probability, a site with f of 0.1 or 0.9 can drift into it; the contrast between the groups is much larger than that noise.
