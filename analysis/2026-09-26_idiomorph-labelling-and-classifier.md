# Mucoromycota idiomorph labelling and the HMM classifier
Status: decided

## Question
Why were 33 sexM-only Rhizopus calls labelled Plus, and can profile HMMs type
sexM vs sexP reliably?

## Data and code version
Classifier commits 7277e10, 67a122c, 6e2f58d, 1a00b0a; margin floor 533e266
(25 bits); fragment scoring 8d80bed.

## Method
1. HMM discrimination test on modelled proteins, leave-one-genus-out
   (`results/2026-09-26_sexMP_hmm/`).
2. Shipped classifier: training = curated refs + sexP-clade members; Zygo 23
   held out (`db/Mucoromycota/classifiers/MAT/manifest.yaml`).
3. Re-scans of 293 Mucoromycota genomes after each change.

## Results
- Cause of the Plus skew: btbA (marked Plus-only) outvoted sexM by bitscore.
  Superseded reading: "sampling bias" (withdrawn 2026-09-26).
- Test: 38/38 correct leave-one-genus-out.
- Classifier (`results/2026-09-26_hmm_classifier/NOTE.md`): LOO 85/85; Zygo 23
  held out 23/23 on scaffolds and contigs.
- Re-score (`results/2026-09-26_mucoro_rescore/NOTE.md`): 34 Plus->Minus,
  13 Minus->Plus, 0 lost; all 33 Rhizopus sexM-only calls become Minus.
- min_margin raised from 11.2 to 25 bits on the FastTree split.

## What changed in detection
08a2616 (btbA non-voting), 6c8050e/4a24ffb (model both genes), classifier
commits above, 4457c3a (btbA in both idiomorphs), 8d80bed (fragments).

## Limits
Validation F3 (`results/2026-09-28_validation_f3_f4/NOTE.md`): typing 0/540
wrong, but the margin does not separate MAT genes from HMG paralogs. See the
MAT-gene gate report.

## Files
`results/2026-09-26_sexMP_hmm/`, `results/2026-09-26_hmm_classifier/`,
`results/2026-09-26_mucoro_rescore/`, `results/2026-09-27_btbA_homothallic/`.
