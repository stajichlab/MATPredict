# MAT-gene gate threshold set by each classifier build
Status: decided

## Question
Should the MAT-gene gate's absolute score threshold come from each classifier
build instead of a fixed 100 bits? A rebuild with new records lowered the
Circinella minor sexP score from 112.5 to 85.9 bits
(results/2026-09-29_circinella_trace/NOTE.md), so a fixed number moves calls
across the gate silently. Curator ruling 2026-09-29: yes.

## Data and code version
- Part A code: commit 9a458de on polish-scope-cuts (PR #9); measured from the
  frozen worktree .claude/worktrees/run-3f755db (the same change before a
  rebase onto aa06975, which only added the two_idiomorphs statement).
- Part B code: curation-umbelopsis at 6e5034b (frozen worktree
  run-6e5034b), rebased onto 9a458de.
- Inputs: 293 BFD Mucoromycota genomes
  (results/2026-09-26_early_diverging/lists/Mucoromycota.tsv, 1 suppressed);
  Zygo 23 on scaffolds and contigs.

## Method
1. Negative set: the 189 non-locus HMG proteins outside the sexM and sexP
   clades ("other_HMG") of the 2026-09-27 FastTree, the same negatives the F3
   validation used (results/2026-09-28_validation_f3_f4/). Written by
   scripts/make_paralog_negatives.py to
   db/Mucoromycota/classifiers/MAT/paralog_negatives.faa. Never used to train.
2. Rule (classifier_build.gate_threshold): the nearest-rank 95th percentile of
   the negatives' best classifier scores against the build's own HMMs. In the
   validation build this percentile was 99.9 bits, the 100 the curator
   adopted (>=100 bits kept 96/108 true full proteins and 9/189 paralogs;
   results/2026-09-28_validation_f3_f4/f3_absolute.txt). Anchoring on the
   negatives makes the threshold move with a rebuild's score calibration.
3. build() writes it to manifest.yaml as `mat_gene_gate` (with the negative
   set's sha256 and how many negatives and held-out positives reach it).
   `--gate-only` adds the block to an existing build without rebuilding.
4. mat_gene_gate.py reads the manifest value; the roster
   `mat_gene_min_score`, then the 100-bit default, are fallbacks only.
   Withheld loci record `mat_gene_min_score_source`.

## Results
| Build | Threshold | Negatives at/above | Held-out positives at/above | LOO |
|---|---|---|---|---|
| PR #9 (076afe4 HMMs) | 99.9 | 10/189 | 75/85 | 85/85 |
| curation-umbelopsis | 98.9 | 9/189 | 76/88 | 88/88 |

Source: manifest.yaml `mat_gene_gate` in each tree.

Part A, 293 genomes vs results/2026-09-28_next_fixes/Mucoromycota_a9fb69c
(compare_vs_a9fb69c.txt): genomes called 253 -> 253; 264 calls unchanged,
0 gained, 0 lost, 0 changed; 18 -> 18 loci withheld by the gate, all with
`mat_gene_min_score_source: manifest` (99.9). Zygo 23: 23/23 locus and
idiomorph on scaffolds and contigs (zygo23_3f755db/score.txt).

Part B, 293 genomes vs results/2026-09-28_umbelopsis_rebased/
Mucoromycota_umb-noguard (compare_umb_vs_noguard.txt): genomes called
253 -> 253; 264 unchanged, 0 gained, 0 lost, 0 changed; 38 -> 38 gate
withholds, all from the manifest (98.9). Zygo 23: 23/23 on both inputs
(zygo23_umb_6e5034b/score.txt).
- Circinella minor GCA_016758965.1: still withheld. Its sexP+rnhA locus
  (JAEPRB010000020.1) scores 85.9 bits, below this build's 98.9, with one
  flank (rnhA). The threshold now tracks the build, but 85.9 is below the
  paralog 95th percentile, so withholding is consistent with the rule and
  with the Lichtheimiaceae ruling (a).
- The four thin Umbelopsis Minus calls are unchanged: U. sp. WA50703 x2
  (margin 37.2), U. vinacea gzUmbVina2 (36.5), U. nana (31.9); all high,
  model-typed. (U. vinacea WA0000051536, the record strain, is 99.5.)

## What changed in detection
- 9a458de (PR #9): classifier_build.gate_threshold, compute_gate,
  update_gate; build() writes mat_gene_gate; scripts/build_idiomorph_hmms.py
  --gate-only; scripts/make_paralog_negatives.py; mat_gene_gate.gate_min_score;
  PR #9 manifest gains mat_gene_gate (99.9); paralog_negatives.faa;
  tests/detect/test_mat_gene_gate.py (+5 tests).
- 6e5034b (curation-umbelopsis): manifest gains mat_gene_gate (98.9).

## Limits
- Classifier rebuilds are not bit-reproducible: rebuilding the PR #9 build
  from the identical training set and tool versions moved held-out
  own-scores by up to 8.0 bits (mean 2.1). Both manifests were therefore
  updated with `--gate-only` on the existing HMMs, not rebuilt. The
  per-build threshold absorbs this drift only if the negatives drift with
  the positives; not tested across many rebuilds.
- The negative set is one set of 189 Mucorales/Umbelopsidales HMG proteins
  from one tree; its composition sets the percentile.
- The 95th percentile was chosen to reproduce the curator's 100 bits on the
  validation build; it is not independently optimised.
- The rule changed no call in either scan, so this measurement shows
  stability, not improved accuracy.

## Curator decisions
- Made 2026-09-29: the threshold comes from each classifier build.
- Open: whether to rebuild classifiers fully when records are added (drift up
  to 8 bits) or keep --gate-only updates between planned rebuilds.

## Files
- results/2026-09-29_gate_threshold/: NOTE.md, jobs.txt, compare.py,
  compare_vs_a9fb69c.txt, compare_umb_vs_noguard.txt,
  Mucoromycota_3f755db/, zygo23_3f755db/, Mucoromycota_umb_6e5034b/,
  zygo23_umb_6e5034b/
