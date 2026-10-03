# Polish-cap V3 and the pre-sign-off regression check (2026-09-30)

## Question
Should the per-family polish cap protect clusters whose unmodelled fragment
already scores strongly on the classifier? And how should any change be checked
before sign-off?

## Data and code version
- Cap replay: frozen `.claude/worktrees/run-capprot` (a2fe1b4 plus an
  uncommitted, env-switched patch). Implementation: branch cap-v3 `cf5cb1e`,
  merged into PR #9 as `52b3ff9`.
- Regression check: PR #9 `8fa2cec` (tool), `ad3ea9b` (validation results).

## Method
- Variants: V1 polish every strong cluster beyond the cap; V2 at most one extra
  per family; V3 rank strong clusters first inside the cap (never adds work).
  Strong = fragment score >= the build's gate and not claimed by a paralog class.
- Stress test: caps 3 and 2 against a cap-off reference on Mucoromycota 293, LCG
  621, Jena 64; Mortierellomycota/Kickxellomycota (290) for false promotions.
- Regression check (`scripts/regression_check.py`, `scripts/run_regression_panel.sh`,
  `testset/regression_panel.tsv`): diff of every locus between baseline and
  candidate runs.

## Results
- Cap 6 (first replay): all variants change only one call, M. pusillus NRRL
  A-13674 (Plus, medium, margin 90.5). 954 of 965 strong clusters were already
  inside the top 6.
- Stress test: of 20 cap-off calls lost by the plain low cap, V1 and V3
  recovered 19, V2 18. Wrong promotions: V3 0, V1 0, V2 1 (M. indicus NRRL
  13082, wrong scaffold). V3 at cap 2 dropped two weak M. genevensis Minus calls
  (margins ~31). Zygo 23 stayed 23/23 at every cap.
- Mortierellomycota/Kickxellomycota: no calls at caps 6 or 3, with or without V3.
- V3 regression check (first real use; baseline ad3ea9b vs cf5cb1e): Ascomycota
  30, Basidiomycota 33, record sources 20, Zygo 23: no changes; Mucoromycota 79:
  12 withheld-only changes; LCG 621: 1 call change (M. pusillus recovered).
- Regression check validation (8eff88e vs a2fe1b4) surfaced all four known cases:
  Circinella minor model change, M. pusillus cap skip, U. nana gain,
  M. griseocyanus downgrade.
- Panel: 166 genomes plus Zygo 23 both inputs (three record source genomes added
  2026-10-01: Schizophyllum GCF_000143185.2, Coprinopsis GCA_016772295.1,
  fen-198 GCA_964255735.1).

## What changed in detection
V3 cap ranking (52b3ff9). Regression tooling (8fa2cec, 806c3c6).

## Limits
- V3 benefit at cap 6 rests on one genome; the stress test adds 19 cases.
- Runtime cost could not be measured reliably (node load).

## Curator decisions
- Made: V3 adopted after the stress test; signed off on the regression summary
  (2026-10-01); regression check required before any record, classifier,
  paralog class, scope or rule change is signed off; three genomes added to the panel.

## Files
`results/2026-09-30_cap_protection/NOTE.md`, `results/2026-09-30_cap_v3/NOTE.md`
(regression/diff/summary.md, lcg_diff/regression_summary.md),
`results/2026-09-30_regression_check_validation/` (PR #9), `docs/regression-check.md`.
