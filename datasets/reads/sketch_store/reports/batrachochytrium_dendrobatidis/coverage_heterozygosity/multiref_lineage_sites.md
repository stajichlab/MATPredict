# Lineage-diagnostic sites from five Bd assemblies (2026-11-15)

Question: the 172,706 sites where JEL423 and CLFT044 differ may hold changes private to those two assemblies. Does the result for the
nine re-mapped isolates (see [site_genotypes.md](site_genotypes.md)) hold on sites that separate the two lineages in every assembly?

Assemblies (user-supplied lineage labels; all from the BFD collection): Brazil lineage CLFT044 (reference), CLFT067, CLFT071; other
lineage JEL423, JAM81 (GCF_000203795.2), RTP6 (GCA_003595275.1). Lineage labels were not checked by this analysis.

## Method
- minimap2 2.30, `-x asm10 -c --cs=short`, each assembly (query) against CLFT044 (target); job 32102347. Same preset as the JEL423 alignment (job 32089404).
- Per site: `jel` (assembly carries the JEL423 base), `clft` (CLFT044 base), `third`, or `na` (no block, or two blocks, cover the site).
- A site is diagnostic when CLFT067 and CLFT071 never carry a non-CLFT044 base, JEL423, JAM81 and RTP6 never carry a non-JEL423 base,
  and each group has at least one assembly that supports it. Script `scripts/bd_ref_lineage_sites.py` (4 tests).
- Read counts at diagnostic sites: the existing depth-controlled counts (12x, job array 32093759) restricted with `bd_site_summary.py --sites`.

## Result 1: most JEL423-vs-CLFT044 sites are not lineage-diagnostic
[tables/multiref_site_states_summary.tsv](tables/multiref_site_states_summary.tsv) (source: `results/2026-11-15_bd_multiref/sites/summary.tsv`, job 32102347)

| assembly | jel | clft | third | na |
|---|---|---|---|---|
| JEL423 | 163,865 | 0 | 0 | 8,841 |
| CLFT067 | 68,558 | 69,441 | 765 | 33,942 |
| CLFT071 | 23,615 | 119,831 | 493 | 28,767 |
| JAM81 | 117,169 | 28,101 | 445 | 26,991 |
| RTP6 | 123,057 | 24,169 | 417 | 25,063 |

Diagnostic sites: 52,555 of 172,706 (30.4%). CLFT067 carries the JEL423 base at 68,558 sites (39.7%), against 23,615 (13.7%) for CLFT071.
So CLFT067 does not carry the CLFT044 base at about 40% of the sites. I did not test why (shared ancestry, admixture or assembly
differences). The Brazil label of CLFT067 was not checked.

## Result 2: the HET-rich result holds on diagnostic sites
[tables/site_genotype_diagnostic.tsv](tables/site_genotype_diagnostic.tsv) (source: `results/2026-11-15_bd_multiref/site_genotype_diagnostic.tsv`)

| isolate | class | sites | clft only | het | jel only | (all 172,706 sites: clft / het / jel) |
|---|---|---|---|---|---|---|
| CLFT044 | CLFT-type | 44,046 | 0.976 | 0.022 | 0.002 | 0.890 / 0.105 / 0.005 |
| CLFT165 | CLFT-type | 38,853 | 0.699 | 0.271 | 0.030 | 0.670 / 0.244 / 0.086 |
| CLFT160 | HET-rich | 40,418 | 0.350 | 0.561 | 0.089 | 0.414 / 0.468 / 0.118 |
| CLFT039 | HET-rich | 40,888 | 0.332 | 0.589 | 0.079 | 0.408 / 0.440 / 0.152 |
| five JEL-type (range) | JEL-type | 32,664 to 34,969 | 0.113 to 0.145 | 0.113 to 0.151 | 0.704 to 0.765 | 0.19-0.21 / 0.15-0.19 / 0.60-0.66 |

CLFT160 and CLFT039 keep a het share (0.56, 0.59) far above every control (0.02 to 0.27) on diagnostic sites. The separation is larger
than on all sites. The JEL-type isolates still show 0.11 to 0.15 CLFT-only and het shares; the cause was not tested (candidates:
mapping errors, residual within-lineage variation, five assemblies being too few to remove every shared site).
Shares are relative and are not mixture proportions.

## Limits
- Lineage membership of the five assemblies comes from the user. Two Brazil assemblies and three other-lineage assemblies are few.
- Sites with `na` in one assembly are not vetoed by it.
- This does not give a third-lineage reference. It does show that the CLFT160/CLFT039 result is not caused by sites private to JEL423 and CLFT044.

## CLFT067 examined (added 2026-11-15, at the curator's request)
Why does the CLFT067 assembly carry the JEL423 base at 39.7% of the divergent sites?
1. Pattern along the genome: the JEL423-base share of CLFT067 is 0.35 to 0.62 (deciles) in 240 of 243 windows of 100 kb with at least 30 usable sites; only 3 windows are below 0.2 and none above 0.8. CLFT071 has 153 of 246 windows below 0.2. JAM81 and RTP6 have 154 of 248 and 170 of 248 windows above 0.8. So the CLFT067 signal is spread over the whole genome, not in blocks. (Source: `results/2026-11-15_bd_multiref/sites/site_states.tsv.gz`.)
2. Its own reads (SRR6375565, 44.9x, mapped to CLFT044 in `bd_coverage_isolates.slurm`; SNP table with QUAL >= 30): at the 68,558 sites where the assembly has the JEL423 base, reads show the JEL423 allele at 20,561 (30.0%; 17,597 het-like); at the 69,441 sites where the assembly has the CLFT044 base, reads show it at 3,686 (5.3%). For CLFT071 (SRR6375567) the figures are 53.0% (assembly JEL base) and 9.2% (assembly CLFT base).
   So the reads and the CLFT067 assembly agree in direction (5.7 times more JEL alleles where the assembly has one), but at 70% of the assembly's JEL-base sites the reads show no JEL allele. The SNP table has no row for sites without a confident variant, so "no JEL allele" means "none called", not "none present".
3. Cause: not determined. Candidates are an assembly that picks one haplotype at random per site in a heterozygous isolate, calls lost to the QUAL and mapping filters, and a different isolate. This analysis cannot choose between them.
4. Sensitivity: with CLFT067 removed from the Brazil set (CLFT071 only), 89,866 sites are diagnostic. Het share (table `tables/site_genotype_diagnostic_no067.tsv`): CLFT160 0.571, CLFT039 0.576; CLFT044 0.079, CLFT165 0.292; five JEL-type 0.092 to 0.124. The CLFT160 and CLFT039 result does not depend on CLFT067.
Wording used in the report: "heterozygous for Brazil-lineage and other-lineage alleles" (measured het share); not "hybrid" and no mixture proportion.
