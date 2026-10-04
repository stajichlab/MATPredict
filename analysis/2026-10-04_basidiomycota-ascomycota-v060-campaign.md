# Basidiomycota and Ascomycota campaign on v0.6.0 (2026-10-03/04)

Status: decided (runs complete; R. toruloides HD ruled 2026-10-04: keep redHD).

## Question
Can MATPredict v0.6.0 run over every BFD Basidiomycota and Ascomycota genome?
What call rates, run times and failures does it give, and what changed against
the 2026-09-26 Basidiomycota run?

## Ruling
J. Stajich, 2026-10-03: run a regression check, then a pilot per phylum, then
all Basidiomycota, then all Ascomycota in waves, on the `exfab` partition.
2026-10-04: fix the code-26 failure on a branch with a PR and v0.6.1; commit
these results with the Mucoromycotina campaign.

## Data
- Code: frozen worktree `run-7c7ed99`. Its `src/` and `db/` are identical to
  tag `v0.6.0`.
- Genomes: BFD `samples.csv`, every row of the phylum with a taxid and a genome
  file. Routing comes from the taxid (lineage routing).
- Partition `exfab`. Accounts: `stajichlab` (384 CPUs, 1 TB shared with all
  other jobs of the user) and `exfab` (32 CPUs, 256 GB).

## Method
1. Regression: the fixed panel (166 genomes) plus Zygo 23 on v0.6.0, diffed
   against the last checked candidate (2026-10-01 scope-only run).
2. Pilots: 163 Basidiomycota genomes, stratified by order (57 orders, >= 2
   each); 246 Ascomycota genomes, stratified by class (30 classes, 92 orders).
   `random.seed(20261003)`.
3. Full runs: Basidiomycota as 3 jobs of 64 CPUs (genomes < 500 Mb) and 1 job
   of 32 CPUs (58 genomes >= 500 Mb, 4 h per-genome limit). Ascomycota as 9
   waves of about 2,157 genomes, 64 CPUs each, 1 h per-genome limit.

## Results

### Regression
0 changed loci in all five panels. Zygo 23: 23/23 on scaffold and contig
inputs. `results/2026-10-03_regression_v060/diff/summary.md`.

### Basidiomycota (3,275 listed)
- 3,270 reports; 5 on the suppress list; 0 failures; 0 timeouts.
- Wall time: 3 jobs of about 1 h, 1 job of 32 min. Median per genome 101 s
  (< 50 Mb) to 716 s (> 500 Mb).

| Subphylum | Genomes | Called v0.6.0 | Without PR-only genomes | Called 2026-09-26 |
|---|---|---|---|---|
| Agaricomycotina | 2,417 | 2,231 (92.3%) | 2,021 (83.6%) | 2,007 (83.0%) |
| Pucciniomycotina | 486 | 391 (80.5%) | 371 (76.3%) | 64 (13.2%) |
| Ustilaginomycotina | 315 | 308 (97.8%) | 301 (95.6%) | 301 (95.6%) |
| Wallemiomycotina | 51 | 51 (100%) | 51 (100%) | 0 |

- Gains come from the records curated since 2026-09-26: `redHD`/`redPR`
  (Sporidiobolales, about 246 genomes), `rustHD` (107), `wallMAT` (51).
- `Basidiomycota:PR` is called in 1,111 more genomes. 1,360 PR calls carry
  `verification: unverified` (admitted only through the strict-CAAX scan). In
  phylum-fallback orders most calls are PR only: Sebacinales 12/12 PR-only,
  Cantharellales 73/90, Trichosporonales 50/60, Auriculariales 18/19. Do not
  count these as MAT loci without the CAAX caveat.
- Apparent losses are mostly merges of one interval into one call: HD into
  Aalpha 96, HD into redHD 15, Schizophyllum Balpha + Bbeta into one B call 40.
- 26 Sporidiobolales (R. toruloides) genomes lost a generic-HD call (HD1/HD2,
  often with MIP1, on the PR contig, about 200 kb from PR). v0.6.0 calls
  `redHD` on a different contig, without MIP1. Which one is the MAT-linked HD
  pair is open (see open questions).
- 2 genomes lost every call: GCA_030280995.1 (Pucciniales, bLocus) and
  GCA_054785265.1 (Kriegeriales, HD).
- Tables: `results/2026-10-03_basidiomycota_v060/genomes.tsv`, `loci.tsv`,
  `by_order.tsv`, `by_order_nonPR.tsv`, `misses_by_order.tsv`, `size_bins.tsv`.

### Ascomycota (19,415 listed)
- 19,380 reports; 15 on the suppress list; 20 failures. All 20 are Alaninales
  (genetic code 26): exonerate has no table 26 and exits 1, which ended the
  genome. Fixed in v0.6.1 (PR #13). Re-run of the 20 on v0.6.1
  (`results/2026-10-04_alaninales_v061/`): 20/20 reports, 0 failures, 30-35 s
  each, miniprot-only models. 4 called (P. tannophilus GCA_001661245.1: MATsc
  MATa, MTL alpha, MAT MAT1-1; 3 others one MATsc call each), 16 uncalled. All
  route `phylum_fallback`: no Alaninales record exists (reference gap).
  Curator (2026-10-04): miniprot models genes less well than exonerate, so pass
  exonerate the NCBI table as a 64-letter string (accepted by stock exonerate
  2.4.0) instead of skipping it. Re-run on that code (`run-b04a77b`,
  `results/2026-10-04_alaninales_v061_exostring/`): 20/20 reports, **15/20
  called** (miniprot-only 4/20); median 457 s per genome (34 s). Regression vs
  v0.6.0: 0 changed loci (`results/2026-10-04_regression_v061s/`). Family labels
  under phylum fallback are not reliable (MATsc and MTL in one genome).
- 17,314 of 19,380 genomes have >= 1 call (89.3%). 848 CPU-h; median 162 s per
  genome; maximum 3,095 s.

| Class | Genomes | Called | Routing |
|---|---|---|---|
| Sordariomycetes | 6,017 | 93.0% | lineage |
| Eurotiomycetes | 3,221 | 96.5% | lineage |
| Pichiomycetes | 2,967 | 89.9% | lineage 2,367 / fallback 575 |
| Dothideomycetes | 2,731 | 89.4% | lineage |
| Saccharomycetes | 2,648 | 85.1% | lineage 2,281 / fallback 364 |
| Lecanoromycetes | 479 | 86.2% | lineage |
| Leotiomycetes | 435 | 92.2% | lineage |
| Dipodascomycetes | 398 | 28.7% | mostly fallback |
| Pezizomycetes | 172 | 79.7% | lineage |
| Orbiliomycetes | 67 | 9.0% | fallback |
| Lipomycetes | 59 | 94.9% | fallback |
| Taphrinomycetes | 22 | 22.7% | lineage |

- 21,011 loci: Ascomycota:MAT 12,628; MATsc 4,954 (includes HML/HMR silent
  copies); MTL 2,799; MATtub 200; PM 138; MATyl 123; mat1 100; mat2 69.
- Confidence: high 13,740, medium 7,112, low 159. Strict pass 20,695, relaxed
  316. `homothallic_candidate` loci: 357.
- Ascomycota:MAT genomes: MAT1-1 only 6,269; MAT1-2 only 5,692; both 204.
- Tables: `results/2026-10-03_ascomycota_v060/genomes.tsv`, `loci.tsv`,
  `by_class.tsv`.

## Decision
Results kept as the v0.6.0 baseline for both phyla. Code-26 fix shipped as
v0.6.1 (PR #13).

## Open
- R. toruloides HD: old generic-HD locus (PR contig) versus v0.6.0 `redHD`
  locus (other contig). Examples: GCA_000222205.2 (AEVR02000017.1:693,823-734,526
  old; AEVR02000008.1:1,252,021-1,268,991 new), GCA_001542265.1,
  GCA_000258745.1. Curator ruling pending.
- Classes with low call rates under phylum fallback (Dipodascomycetes,
  Orbiliomycetes, Trigonopsidomycetes) are reference gaps, not measured
  absence.
- The genome with 11 Ascomycota loci has not been examined.

## Files
- `results/2026-10-03_regression_v060/`
- `results/2026-10-03_basidiomycota_pilot/`, `results/2026-10-03_ascomycota_pilot/`
- `results/2026-10-03_basidiomycota_v060/`, `results/2026-10-03_ascomycota_v060/`
  (`reports_all.tar.zst` = every `detection_report.yaml` and `wall_seconds`;
  per-genome GFF3 and FASTA are not committed)
- `results/2026-10-04_alaninales_v061/`
