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
2. (Superseded by the read-count test in the last section: the SNP-table figures below count sites with a called SNP, not allele fractions.) Its own reads (SRR6375565, 44.9x, mapped to CLFT044 in `bd_coverage_isolates.slurm`; SNP table with QUAL >= 30): at the 68,558 sites where the assembly has the JEL423 base, reads show the JEL423 allele at 20,561 (30.0%; 17,597 het-like); at the 69,441 sites where the assembly has the CLFT044 base, reads show it at 3,686 (5.3%). For CLFT071 (SRR6375567) the figures are 53.0% (assembly JEL base) and 9.2% (assembly CLFT base).
   So the reads and the CLFT067 assembly agree in direction (5.7 times more JEL alleles where the assembly has one), but at 70% of the assembly's JEL-base sites the reads show no JEL allele. The SNP table has no row for sites without a confident variant, so "no JEL allele" means "none called", not "none present".
3. Cause: not determined. Candidates are an assembly that picks one haplotype at random per site in a heterozygous isolate, calls lost to the QUAL and mapping filters, and a different isolate. This analysis cannot choose between them.
4. Sensitivity: with CLFT067 removed from the Brazil set (CLFT071 only), 89,866 sites are diagnostic. Het share (table `tables/site_genotype_diagnostic_no067.tsv`): CLFT160 0.571, CLFT039 0.576; CLFT044 0.079, CLFT165 0.292; five JEL-type 0.092 to 0.124. The CLFT160 and CLFT039 result does not depend on CLFT067.
Decisions of the curator (J. Stajich, 2026-11-16, answered in the session):
- Primary diagnostic set: the 89,866 sites defined without CLFT067 (`diagnostic_sites_no067.tsv`; table `site_genotype_diagnostic_no067.tsv`). The 52,555-site set (with CLFT067) is the sensitivity check. The CLFT067 assembly disagrees with its own reads (see the allele-fraction section), so it does not define diagnostic sites.
- Wording for CLFT160 and CLFT039: "heterozygous for Brazil-lineage and other-lineage alleles" (measured: 51.0% and 52.6% of diagnostic sites have a JEL allele fraction of 0.3 to 0.7; CLFT044 control 6.6%). A mixed culture of two isolates is not excluded by these data, because it would also give a peak near 0.5. No "hybrid", no mixture proportion, no ploidy claim.

## CLFT067 allele-fraction test (added 2026-11-16; answers "is CLFT067 a Brazil x GPL hybrid?")
Question. The assembly CLFT067 (GCA_036289345.1) carries the JEL423 base at 68,558 of 172,706 divergent sites (39.7%), in 240 of 243 100-kb windows. Do the reads of CLFT067 show the same mix?
Method. SRR6375565 (CLFT067; 10,836,360 pairs, matched ENA; median depth 51.1x) was mapped to CLFT044 (bwa-mem2, duplicates removed, no downsampling; job 32112178, `scripts/bd_site_genotype.slurm`). Reads with the JEL423 or CLFT044 base were counted at each divergent site (`bd_site_allele_counts.py`; q >= 20, Q >= 20; this counts reads, so a missing SNP row is not read as "absent"). The profile (`scripts/bd_allele_fraction_profile.py`) uses the 89,866 diagnostic sites defined without CLFT067 (`diagnostic_sites_no067.tsv`) and sites with at least 8 reads. Controls are the 12x counts of CLFT044, CLFT160 and CLFT039. Tables: `tables/clft067_allele_fraction_hist.tsv`, `tables/clft067_allele_fraction_windows.tsv`.
Results (f = JEL reads / all reads; "mid" = 0.3 <= f <= 0.7).
| Isolate | Sites | f < 0.1 | f 0.3-0.7 (mid) | Mid expected if heterozygous at this depth (binomial model) |
|---|---|---|---|---|
| CLFT067 (51x) | 83,754 | 72,540 (86.6%) | 10.3% | 98.7% |
| CLFT044 (12x, reference lineage) | 65,933 | 59,339 (90.0%) | 6.6% | 86.1% |
| CLFT160 (12x) | 57,555 | 12,070 (21.0%) | 51.0% | 85.1% |
| CLFT039 (12x) | 59,041 | 11,302 (19.1%) | 52.6% | 86.0% |
1. All 172,706 divergent sites, CLFT067 reads (at least 5 reads, 151,114 sites): 83.1% carry the CLFT044 base, 14.8% are mixed (f 0.2 to 0.8), 2.1% carry the JEL423 base.
2. At the diagnostic sites defined without CLFT067, the CLFT067 reads look like the Brazil lineage at most sites. They do not show the genome-wide mid-range peak that CLFT160 and CLFT039 show. In 240 of 240 100-kb windows (at least 20 sites) the median f is 0.000. In CLFT160, 177 of 232 windows have a median f between 0.3 and 0.7.
3. CLFT067 has more mid-range sites than CLFT044 (10.3% against 6.6%; the depths differ, 51x against 12x, so the two values are not exactly comparable). 65 windows have a mid share above 0.2 and none above 0.5 (maximum 0.42). A minor fraction, a few regions or a small mixed component are possible; I did not test which.
What this means.
- The reads of SRR6375565 do not give an F1-like pattern (a peak at f = 0.5 at most sites). This statement is narrowed in the follow-up below: the reads do carry a minority JEL signal at the sites where the assembly has the JEL base. The earlier "30%" figure was a SNP-table detection fraction and is replaced by this test.
- The assembly and the reads disagree: the assembly has the JEL base at 39.7% of sites; the reads have it at 2.1%. The cause is undetermined. Candidates: the assembly includes other-lineage sequence (mis-assembly or contamination), the assembly and SRR6375565 come from different cultures of CLFT067, or the mapping to CLFT044 hides those sites. I did not test them.
- CLFT067 stays in the Brazil set only on the strength of its reads. The assembly alone would not support that.
Limits. One read set per isolate; depths differ between CLFT067 (51x) and the controls (12x); 12x-downsampled CLFT067 was not run.

Follow-up, same day: reads split by the state of the CLFT067 assembly (corrects the summary above).
The diagnostic set above is defined by other assemblies. A second split uses the CLFT067 assembly's own base at each of the 172,706 divergent sites (same reads, at least 8 reads per site).
| CLFT067 assembly base at the site | Sites with reads | Reads with JEL base (pooled) | f < 0.1 | f 0.1-0.9 | f > 0.9 |
|---|---|---|---|---|---|
| JEL423 base | 63,916 | 17.3% | 42,260 | 19,330 (30%) | 2,326 |
| CLFT044 base | 63,992 | 2.9% | 59,719 | 3,950 (6%) | 323 |
- At sites where the assembly has the JEL base, 30% of sites are mixed in the reads and 17.3% of reads carry the JEL base. At sites where the assembly has the CLFT044 base, 6% of sites are mixed. So the reads do contain a minority JEL component that the assembly sits on.
- By 100-kb window (237 windows with at least 20 such sites) the pooled JEL read share at assembly-JEL sites has median 0.19 (quartiles 0.05 and 0.28, range 0.00 to 0.73). 77 windows are below 0.1, 85 are at least 0.25 and 11 are at least 0.4. The highest windows are on JAQIHF010000002.1 and JAQIHF010000001.1.
- This pattern fits a mixed or heterozygous component at a minority fraction in some regions, with the assembly carrying the JEL base. It does not fit a clean F1 (peak at 0.5 everywhere) and it does not fit a pure Brazil genome. A hybrid or recombinant origin, a mixed culture and a haplotype choice in the assembly all remain possible; these data do not separate them.
- The earlier PR comment of 2026-11-16 said the reads "do not show a genome-wide Brazil x GPL mix". That is correct for the diagnostic set and for an F1-like pattern, but incomplete: see this follow-up.

## CLFT067 reads mapped to the CLFT067 assembly (added 2026-11-16, job 32114411; decision 7)
Test: SRR6375565 (ENA read count matched) mapped with bwa-mem2 to its own assembly (contigs JAZDSG01..., 264 windows of 100 kb). Flagstat: 97.77% mapped (`tables/clft067_selfmap_flagstat.txt`). Script `scripts/bd_self_map_windows.slurm`; comparison `scripts/bd_selfmap_vs_jel_share.py` (sites moved from CLFT044 to assembly coordinates by linear interpolation inside primary PAF blocks; approximate).
Results (233 windows with at least 50 diagnostic-state sites; `tables/clft067_selfmap_windows_vs_jel_share.tsv`):
- Median window depth 77.6x (10th to 90th percentile 68.7-151x). This is higher than the 51x quoted earlier, which came from the CLFT044 mapping; I did not resolve the difference.
- Read mismatch (edit distance per aligned base, MAPQ >= 20): genome 0.0076, window median 0.0068, range 0.0029-0.053.
- The assembly's JEL-base share by window has median 0.51 (range 0.13-0.94). Only one window has share > 0.9 and none < 0.1. So the JEL-state sites are spread through all windows, not held in a few JEL-like blocks.
- Spearman: JEL share vs depth 0.20, JEL share vs edit rate 0.18, depth vs edit rate 0.59. The depth-mismatch link points to collapsed or repeat windows; the JEL-share links are weak.
Reading: a window-level mix of two haplotypes would show higher depth and edit rate where the assembly holds the JEL haplotype. Only weak correlations are seen. The test does not separate a hybrid, a mixed culture and an assembly of mixed haplotypes. Open.
Limit: the interpolation lift and the 100 kb window are coarse; the window JEL share of 0.51 here differs from the 39.7% genome share quoted above because it is the median of windows with at least 50 sites.
