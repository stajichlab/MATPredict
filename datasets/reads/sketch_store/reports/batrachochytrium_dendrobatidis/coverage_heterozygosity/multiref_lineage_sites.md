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
