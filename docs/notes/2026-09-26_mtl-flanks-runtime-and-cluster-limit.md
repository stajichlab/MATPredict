# MTL flanks: are they useful, what do they cost, and can polishing be capped?

2026-09-26. Data: `results/2026-09-26_flank_ablation/` (60 Serinales genomes x 7
variants, `summary.txt`), `results/2026-09-26_cluster_limit/` (`study_output.txt`,
`study2_output.txt`), and the 2,368-genome Serinales baseline scan
`results/2026-09-25_serinales_all_8dec9e1/`.

## 1. The flanks do the work; MTLalpha2 does not change a call

Commit `2c6c3a7` added PAP1/OBP1/PIK1 as optional flanks and MTLalpha2 as an
optional core gene in one step. The ablation separates them. Sample: 20 C.
albicans, 17 C. auris, 23 other Serinales across 20 genera, stratified and
seeded.

| variant | C. albicans | C. auris | other Serinales | same genotype as full roster |
|---|---:|---:|---:|---:|
| V0 baseline | 10/20 | 15/17 | 9/23 | 29/60 |
| V1 + MTLalpha2 only | 19/20 | 15/17 | 9/23 | 46/60 |
| V2 + flanks only | 20/20 | 16/17 | 20/23 | **60/60** |
| V3 both (`2c6c3a7`) | 20/20 | 16/17 | 20/23 | 60/60 |
| V4 flanks, miniprot-only, no alpha2 | 20/20 | 16/17 | 20/23 | 60/60 |
| V5 flanks, miniprot-only, + alpha2 (`882aa01`) | 20/20 | 16/17 | 20/23 | 60/60 |

* Flanks alone reproduce the full result. MTLalpha2 alone recovers C.
  albicans (it gives alpha a second core gene) but not the rest of Serinales
  (9/23 vs 20/23) nor the C. auris alpha call.
* V4 vs V5: identical genotypes in 60/60. MTLalpha2 appears in 23 V5 calls as
  supporting evidence. Kept as a reference, per the curator.
* 0 calls in any variant rest on flank models without a modelled core gene.
* Flank-only clusters (4,937 in 108 genomes of the partial scan, ~46 per genome)
  are OBP1/PIK1/PAP1 paralog and background hits with no MTL gene; the evidence
  floor admits none of them. They cost tblastn time, not polish time, and never
  become calls.
* C. albicans genotypes: only 6/20 a+alpha, 10 alpha-only, 4 a-only. Wild
  isolates are mostly a/alpha heterozygous; the likely explanation is that
  these diploid assemblies collapse the two MTL haplotypes. **Unverified**; it
  needs read-level checks.

## 2. Runtime: exonerate on the long flanks was the cost

Profile, C. albicans 3153A: baseline 15 s; with flanks 40 s, of which 32 s was 14
exonerate calls on PAP1/OBP1/PIK1 (~2.3 s each) and 0.5 s miniprot on the same
genes; genome-wide tblastn 3.4 s.

Fix (`882aa01`): a roster gene may declare `polish: miniprot`, which skips
exonerate for that gene; set for the three flanks only. On the 60 genomes, V5
matches the pre-implementation driver (V3m) in 60/60 genotypes, and median
wall per group fell from 150-204 s (V3) to 42-66 s (V5). These timings were
taken on a node shared with other jobs: relative, not absolute.

Rejected, curator ruling: searching flanks only near core-gene hits. It would
change the search database per genome and could miss loci, for a modest gain.

## 3. Capping polish per genome: projected to halve the slow panels

In the slow panels most admitted clusters never become a call, and wall time is
close to linear in the admitted-cluster gene load:

| panel | admitted clusters / genome (median) | became a call | wall ~ a + b x genes | R^2 |
|---|---:|---:|---|---:|
| Dothideomycetes | 26 | 3% | 127 + 8.75 | 0.81 |
| uncurated Pezizomycotina | 19 | 4% | 148 + 9.83 | 0.65 |
| Orbiliomycetes | 18 | 2% | 43 + 8.24 | 0.70 |
| Dipodascomycetes etc. | 32 | 2% | 54 + 6.51 | 0.92 |
| other Saccharomycetes | 27 | 2% | 20 + 8.54 | 0.83 |
| Serinales baseline | 3 | 22% | 11 + 0.72 | 0.40 |

Limit: polish at most N admitted clusters, ranked by what is known BEFORE
polishing (distinct genes, then best identity, then hit count). Ranking PER
FAMILY matters: a phylum-fallback genome is searched against many families,
and a genome-wide top-6 loses 8/80 calls there against 1/80 per family.

Per-family top-6:

| panel | calls lost | projected median wall |
|---|---:|---|
| Dothideomycetes | 0/91 | 622 -> 284 s |
| uncurated Pezizomycotina | 0/88 | 620 -> 315 s |
| Orbiliomycetes | 0/28 | 336 -> 191 s |
| Dipodascomycetes etc. | 0/63 | 551 -> 268 s |
| other Saccharomycetes | 1/80 | 516 -> 284 s |
| Saccharomyces (3 cassettes/genome) | 1/3,101 | no wall recorded |
| Serinales baseline | 3/1,327 | ~14 s already |

**This is a projection from existing runs, not a measured result.** The loss
count assumes nothing else changes when clusters are skipped (the relaxed pass,
for example, keys off the strict result). Not implemented. Next step: implement
behind a parameter and measure on the same panels.
