# Would stricter flank support withhold the paralog-type label errors?
Status: open

## Question
Five labelled LCG strains disagree with their strain-name label in a way that
looks like an HMG paralog call (results/2026-09-29_strain_labels_and_absidia/).
Can a stricter MAT-gene gate withhold them without losing correct calls?

## Data and code version
- Gate logic: src/MATPredict/detect/mat_gene_gate.py at remote
  polish-scope-cuts 076afe4, reproduced from report fields in replay.py.
- Reports (read-only): LCG results/2026-09-28_lcg_holdout/runs2 (076afe4);
  Mucoromycota 293 results/2026-09-28_next_fixes/Mucoromycota_a9fb69c;
  Jena results/2026-09-28_mucor_jena_holdout/runs_scaffolds;
  Zygo 23 results/2026-09-28_next_fixes/zygo23_a9fb69c/scaffold.
- Per-build gate threshold: gate-threshold branch commit 3f755db writes
  `mat_gene_gate.min_score: 99.9` for the PR #9 classifier
  (.claude/worktrees/gate-threshold/db/Mucoromycota/classifiers/MAT/manifest.yaml).
- Labels: results/2026-09-29_strain_labels_and_absidia/lcg_label_scores.tsv
  (clean tier; the 3 disputed labels excluded).

## Method
1. Reproduce the current gate route (absolute score, flank support, fail) of
   every reported classified call. No reported call falls in "fail", so the
   replay reproduces the gate (replay_output.txt, top).
2. Variants:
   - V1 (= V3): flank route needs >=2 flanks at >=40%, at least one not rnhA.
   - V2(X): partial_locus below the absolute threshold withheld unless
     margin >= X (X = 75, 100).
   - V2b(X), added: any partial_locus withheld unless margin >= X.
3. Count withheld calls per set, targets withheld, Zygo, label agreement.

## Results
| Variant | Targets withheld (of 5) | Withheld LCG / Mucoro293 / Jena | Zygo 23 | Clean labels agree/disagree/uncalled |
|---|---|---|---|---|
| current | – | – | 23/23 | 10 / 5 / 4 |
| V1 = V3 | 0 | 0 / 0 / 0 | 23/23 | 10 / 5 / 4 |
| V2 X=75 | 0 | 12 / 6 / 1 | 23/23 | 10 / 5 / 4 |
| V2 X=100 | 0 | 12 / 9 / 1 | 23/23 | 10 / 5 / 4 |
| V2b X=75 | 3 | 46 / 11 / 2 | 23/23 | 8 / 2 / 9 |
| V2b X=100 | 3 | 56 / 15 / 4 | 23/23 | 8 / 2 / 9 |

- Why V1 and V2 miss the targets: the three partial targets (Circinella
  angarensis RSA_198_Plus 132.5 bits, C. umbellata RSA_505_Plus 134.4,
  Thamnostylum repens RSA_459_Plus 109.3) pass the gate on the ABSOLUTE route
  (model score >= 100), not on flank support. Backusella NRRL_6044_Plus
  (mat_locus, 155.2 bits, glrA+algA) and Pilaira RSA_1997_Plus (mat_locus,
  flank route with algA, tptA, rnhA) pass every variant.
- The per-build threshold (99.9) changes nothing for these runs.
- V2 withholds only fragment-typed low calls (mostly Mucor genevensis x7
  Minus, reported homothallic; Umbelopsis fragment calls; gzMucHiem1).
- V2b withholds the 3 partial targets but also 2 label-agreeing calls
  (C. angarensis RSA_618- and Fennellomyces RSA_1418-, both Minus) and
  ~40 other partial calls (13 Lichtheimiaceae/Circinella in LCG).
- Pattern (inferred, not tested): the partial sexM+rnhA calls carry identical
  classifier scores across strains labelled OPPOSITE mating types and across
  genera — C. umbellata NRRL 1366 and RSA_505_Plus both 134.4/69.0 (also
  Rhizopus arrhizus NRRL 1470); C. angarensis RSA_618- 133.5 vs RSA_198_Plus
  132.5; Helicostylum and Thamnostylum 124.3/63.0; three Syncephalastrum
  102.6/28.0. A sequence shared by Plus- and Minus-labelled strains is not
  idiomorph-specific; this looks like a conserved sexM-like HMG gene next to
  rnhA present in both mating types.

## What changed in detection
Nothing. Read-only replay.

## Limits
Label truth is 18 clean LCG strains. V2b's collateral losses are unverified
beyond the 2 labelled strains. Identical scores suggest a shared gene but no
alignment was done.

## Curator decisions
Open: none of the flank variants reach the targets. The candidate next test is
whether the called "sexM" at these partial loci is identical between opposite-
labelled strains (align the models) — a shared-in-both-mating-types check.

## Files
replay.py, replay_output.txt (this folder).
