# Explicit MAFFT L-INS-i and the candidate deterministic rebuild

Status: decided (code: explicit L-INS-i, no trimming); open (curator approval
to ship the candidate full rebuild).

## Question
1. Set the classifier aligner to MAFFT L-INS-i explicitly instead of `--auto`.
2. With that setting, what would a deterministic full rebuild of the
   Mucoromycota sexM/sexP classifier change on Mucoromycota, Zygo 23, the LCG
   and Jena held-out sets?

## Data and code version
- Code: PR #9 `polish-scope-cuts` at f084d2f (this change) on top of 64971b6.
- Candidate: built with scripts/build_idiomorph_hmms.py (full build) into
  candidate/ (not db/). Measured from frozen worktree `run-algcand`
  (branch aligner-candidate-measure, commit 4f29df3 = f084d2f + candidate HMMs,
  measure-only).
- Baseline: the shipped classifier, PR #9 R4 runs (41bd471),
  results/2026-09-29_r4_paralog/.
- Inputs: Mucoromycota 293 (results/2026-09-26_early_diverging/lists, taxid
  routing, same as baseline); Zygo 23 on scaffolds and contigs; LCG
  Mucoromycotina 621 and Jena 64 scaffolds with `--phylum Mucoromycota`.
- SLURM jobs in jobs.txt.

## Method
1. MAFFT_OPTIONS = `--localpair --maxiterate 1000 --quiet --thread 1`; test
   that `--auto` is absent (tests/detect/test_classifier_build_deterministic.py).
2. Full build of Mucoromycota:MAT (sexP, sexM, P1) into candidate/.
3. Checksums against the earlier `--auto` candidate
   (results/2026-09-29_deterministic_build/candidate/) and the shipped HMMs.
4. The same five runs as R4 with only the classifier files swapped; compare.py.

## Results
- Tests: 938 pass.
- Candidate build: LOO 85/85, worst correct margin 21.8 bits; gate threshold
  100.3 bits (shipped 99.9).
- Checksums: the explicit L-INS-i candidate is byte-identical to the earlier
  `--auto` candidate (sexP 98f042b0..., sexM a8d08271..., P1 a802fe1a...),
  confirming `--auto` had chosen L-INS-i. Both differ from the shipped HMMs.

| Set | Genomes | Called (shipped -> candidate) | Lost | Gained | Changed |
|---|---|---|---|---|---|
| Mucoromycota | 293 | 253 -> 253 | 0 | 0 | 0 |
| LCG Mucoromycotina | 621 | 536 -> 536 | 0 | 0 | 0 |
| Jena scaffolds | 64 | 61 -> 61 | 0 | 0 | 0 |

- "Changed" covers any change of idiomorph or confidence at the same locus.
- Withholdings identical: paralog_class 13 / 33 / 8 and mat_gene_gate
  7 / 22 / 1 (Mucoromycota / LCG / Jena) in both arms.
- Classifier margins on 893 matched calls moved by -8.1 to +11.3 bits
  (mean |d| 2.62); no call crossed the 25-bit typing floor or the gate.
- Zygo 23: 23/23 locus and idiomorph on scaffolds and on contigs.
- LCG clean strain labels (n=17; disputed, misidentified and leaked strains
  excluded): agree 10, disagree 3, uncalled 4 in both arms.

## What changed in detection
- f084d2f: explicit L-INS-i in src/MATPredict/detect/classifier_build.py;
  shipped HMMs unchanged.

## Limits
- Mucoromycota classifier only.
- The margin shifts (up to 11.3 bits) are the one-time replacement of the
  non-deterministic shipped build; later rebuilds of the same inputs are
  byte-identical.

## Curator decision
Open: ship the candidate as db/Mucoromycota/classifiers/MAT/ (0 call changes
across 978 genomes and Zygo 23), or keep the shipped HMMs.

## Files
compare.py, compare_output.txt, changes.tsv (empty of changes), jobs.txt,
candidate/, run_lcg_chunk.slurm, run_jena.slurm, Mucoromycota_4f29df3/,
lcg_runs/, jena/, zygo23_4f29df3/.
