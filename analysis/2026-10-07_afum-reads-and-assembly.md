# Aspergillus fumigatus: MAT typing from 304 assemblies and 331 read sets
Status: open (first full profile; decisions below are open)

## Question
How well are MAT1-1 and MAT1-2 retrieved from assemblies (`detect`) and from raw reads (`matpredict reads-type`) in a large *A. fumigatus* population set, and how do the two compare? Tool and method: [2026-10-06_fola-reads-type.md](2026-10-06_fola-reads-type.md) (Fola benchmark, mixes, blastx and recruit tests).

## Data and code version
- Code: branch `afum-reads-type` (stacked on `reads-type`, PR #45), `matpredict reads-type` with `--min-unique-run 100`. `detect`: frozen worktree `.claude/worktrees/run-2d1860a` (main at PR #41), `--taxid 746128`.
- Database: the same tree as the `detect` code. The panel uses records `746128_a1163_MAT_MAT1-1` and `746128_af293_MAT_MAT1-2`.
- Inputs (`results/2026-10-07_afum_reads/strain_map.tsv`):
  - Assemblies: `/bigdata/stajichlab/shared/projects/Afumigatus_pangenome/scaffolded/genomes/*.sorted.fasta`, 304 files.
  - Reads: `/bigdata/stajichlab/shared/projects/Afum_popgenome/asm/input/`, 331 read sets from `asm/samples.dat` (column 1 = read-file prefix, column 2 = strain = genome name). 311 are paired and 20 are single-end (`_R1` only; typed as one file).
  - 303 genomes have reads (Af10 has none); 28 read sets have no genome. 20 DMC and DMC2 pairs (for example DMC_AF100-1_3 and DMC2_AF100-1_3) give identical read results, so they are duplicates under two names.

## Method
1. Panel (`panels/aspergillus_fumigatus/`): the curated database records `746128_a1163_MAT_MAT1-1` (A1163, DS499596.1, 13,260 bp) and `746128_af293_MAT_MAT1-2` (Af293, NC_007196.1, 13,096 bp), whole locus segments with the flank genes. The two segments align at 99.3% over the flanks. The idiomorph-specific regions are MAT1-1 positions 4,576 to 6,664 (2,089 bp) and MAT1-2 positions 6,225 to 8,525 (2,301 bp) (`specific_regions.fasta`).
2. Assembly truth (`scripts/assembly_idiomorph_blast.py`, `blast/assembly_truth.tsv`): BLASTN of the two specific regions against each assembly (identity >= 90%, e-value 1e-10); present if >= 0.80 of the region is covered, absent if <= 0.20, otherwise partial. Result for 304 genomes: MAT1-1 only 160, MAT1-2 only 128, both 8, partial 7 (excluded from agreement counts), neither 1 (Moose_3-1.1).
3. `detect` (job 29605756) on all 304 assemblies; `reads-type` (jobs 29605757 first run, 29606456 with the filter) on 331 read sets, first 8M reads each, panel `panels/aspergillus_fumigatus/`; `scripts/afum_compare.py` joins the three calls per strain.

## Results
**A panel flaw found and fixed.** With the first run (all unique k-mers), reads agreed with the assembly truth in 288 of 296 strains (97.3%), and 6 MAT1-2 assemblies were called `both`. In MAT1-2 strains the MAT1-1 "unique" k-mers had median breadth 0.351 and about half the MAT1-2 depth, as high as the signal in those 6 strains. Cause: the A1163 and Af293 flanks differ by SNPs, and every k-mer that spans a SNP counted as idiomorph-unique, and matches any strain with that flank allele. Fix (`Panel.from_sequences(min_run=100)`, option `--min-unique-run`, default 100): keep a unique k-mer only if it lies in a run of at least 100 consecutive unique positions (a SNP gives a run of at most 31). It removes 1,175 of 3,739 MAT1-1 and 1,168 of 3,575 MAT1-2 unique k-mers in this panel (31%), and 111 of 4,828 and 92 of 4,711 in the Fola panel (2%). Tests: 3 new, 47 pass with the smoke test. The rule was chosen after seeing the first result, so the 99.0% below is not an independent test of it. Fola is unchanged: 147/148 (`v4_concordance.tsv`).

Results with the filter (`comparison.tsv`; excluding the 7 partial and the 28 reads-only strains):
| Method | Agree with assembly truth |
|---|---|
| reads-type, 8M reads | 293/296 (99.0%); without the 20 DMC2 duplicates 280/283 |
| `detect` on assemblies (frozen run-2d1860a) | 288/297 (97.0%); without the DMC2 duplicates 275/284 |
- reads-type: MAT1-1 160/160, MAT1-2 126/127, both 7/8. Differences: AF100-12_2 (below), F18149-Manchester (below), Moose_3-1.1 (`low_depth`).
- `detect`: all 288 assemblies with one idiomorph got the right idiomorph (160 MAT1-1, 128 MAT1-2), 0 missed. All 8 assemblies with both regions got a single call (5 MAT1-1, 3 MAT1-2); the second idiomorph is not reported. All 304 genomes had at least one locus; 185 loci high and 120 medium confidence; no run failed.
- All reads (331): MAT1-1 174, MAT1-2 135, both 18, low_depth 4. Reads-only strains (28): MAT1-1 11, MAT1-2 7, both 7, low_depth 3.

The 8 assemblies with both regions (read depths, MAT1-1 / MAT1-2, k-mer units): 08-36-03-25 4.3/14.0, AF100-12_5 19.3/3.8, Afu_343_P_11 15.4/3.3, B7586_CDC-30 26.9/4.2, DMC_AF100-1_18 22.2/3.8, IFM_59359 24.0/16.4, IFM_61407 15.8/14.0: all `both` from reads. F18149-Manchester: reads call MAT1-1 (16.5) with MAT1-2 at 1.4 (breadth 0.71, flagged `trace_MAT1-2`), so the reads support the second region only weakly although the assembly holds both. Two of the three Lofgren et al. putative hybrids (IFM_59359, IFM_61407) are `both` from reads and from the assembly regions; `detect` reports MAT1-2 only (medium) for both. The third (AF100-1_3) has no assembly in this folder; its reads (DMC_AF100-1_3, DMC2_AF100-1_3) are `both`, at 13.3 / 13.4.

Reads show a second idiomorph that the assembly lacks: AF100-12_2 (assembly MAT1-2 only; reads MAT1-1 depth 3.5 against MAT1-2 18.8, MAT1-1 breadth 0.546). It is marginal: in the 126 MAT1-2 strains called MAT1-2 the MAT1-1 breadth has median 0.144 and maximum 0.443, and all 126 carry a `trace_MAT1-1` flag. In the 160 MAT1-1 strains the MAT1-2 breadth has median 0.000, maximum 0.514, and 46 carry a flag.

`detect` findings that are not read-typing results:
- 168 of the 169 MAT1-1 loci are labelled `homothallic_candidate`, and none of the 135 MAT1-2 loci are (`mat_locus`). The A1163 database record carries a degenerate MAT1-2-1 remnant as a core gene, and the class rule (`idiomorph.py`, both idiomorphs' core genes on one locus) fires on it. For *A. fumigatus* this label is therefore not evidence of homothallism. I did not change the record or the rule.
- Moose_3-1.1: `detect` reports MAT1-2 (high, scaffold_3:476,818-485,874) and MAT1-1 (medium, `partial_locus`, scaffold_3:1,319,318-1,326,086), 840 kb apart. Neither specific region matches the assembly at 90% identity, and its reads have shared depth 0.02, so the reads are not from an *A. fumigatus*-like genome. The species of this sample was not checked.

## What changed in detection
- New in `reads-type`: `--min-unique-run` (`src/MATPredict/reads/panel.py`). No change to `detect`, the database, or any classifier.

## Limits
Limits: the assembly truth uses regions from the same two reference strains as the panel; it is independent of the reads but not of the reference choice. 20 duplicate read-set pairs inflate n. Only the first 8M reads were used; single-end samples (20) have half the reads of a pair at the same cap. A `both` call does not separate a heterokaryon, a diploid, a mixed culture or contamination. Reads-only strains have no truth. Cut-offs (breadth 0.50, relative depth 0.10, `min_run` 100) are untested outside Fola and this set.

## Curator decisions
Made: none yet.
- Open: keep `--min-unique-run` at 100.
- Open: The `homothallic_candidate` label on every *A. fumigatus* MAT1-1 locus: change the record (drop the remnant from the core genes), the rule, or leave and note.
- Open: Species of Moose_3-1.1.
- Open: Whether the second idiomorph in AF100-12_2, AF100-12_5, F18149-Manchester is a mixed sample (the AF100-12 series suggests one lab series) and who holds the cultures.
- Open: Whether `detect` should get the report-only "unlinked second idiomorph" check noted in the spec.

## Files
- `results/2026-10-07_afum_reads/`: `strain_map.tsv`, `assembly_truth.tsv`, `comparison.tsv` (per strain: truth, `detect` call, reads call, depths), `reads_type_v2.tsv`, `detect.slurm`, `reads_type_run_v2.slurm`, `blast_run.slurm`.
- `panels/aspergillus_fumigatus/` (MAT1-1.fasta, MAT1-2.fasta, specific_regions.fasta).
- `scripts/assembly_idiomorph_blast.py`, `scripts/afum_compare.py`.
- Full `detect` reports (304 genomes) stay in the main checkout, `results/2026-10-07_afum_reads/detect/runs/` (not tracked).
