# Ship the deterministic explicit-L-INS-i Mucoromycota classifier

## Question
Adopt the deterministic candidate build as the shipped classifier in
db/Mucoromycota/classifiers/MAT/ (curator ruling 2026-09-29)?

## Data and code version
- Code: polish-scope-cuts f084d2f (explicit MAFFT L-INS-i). Worktree
  .claude/worktrees/ship-classifier.
- Training: curated Mucoromycota records + training_extra.faa (unchanged);
  training_sha256 ccfbf37a1ed52dc3d21e864689f33eb233cd60f76ea3917c3899bd85ad055340.

## Method
Full build: `scripts/build_idiomorph_hmms.py --db-root db --family
Mucoromycota:MAT --out db/Mucoromycota/classifiers/MAT`. Checksums against
results/2026-09-29_aligner_default/candidate/.

## Results
- sexM, sexP and paralogs/P1.hmm are byte-identical to the candidate
  (sha256 prefixes a8d0827177fed60f, 98f042b0dfb10269, a802fe1a9719ae31).
- Manifest identical to the candidate apart from notes/build date. It records
  pyhmmer 0.12.3, HMMER 3.4, MAFFT 7.526, options `--localpair --maxiterate
  1000 --quiet --thread 1`, seed 42, training checksum, gate.
- LOO 85/85, worst correct margin 21.8 bits; recommended min_margin 10.9
  (roster min_margin stays 25). Gate 100.3 bits: 9/189 paralog negatives and
  75/85 held-out MAT proteins at or above.
- Replay of this exact classifier (results/2026-09-29_aligner_default/): 0
  calls lost/gained/changed on Mucoromycota 293, LCG 621, Jena 64; Zygo 23
  23/23 locus and idiomorph on both inputs.
- Tests: 938 passed.

## What changed in detection
Commits 9c39e2c (classifier) and 8eff88e (results) on polish-scope-cuts,
pushed fast-forward.

## Limits
Replay done on the candidate files before shipping; the shipped files are
byte-identical, so no separate rerun was made.

## Curator decisions
Made: adopt the candidate (2026-09-29). Open: none for this step.

## Files
db/Mucoromycota/classifiers/MAT/{manifest.yaml,sexM.hmm,sexP.hmm,paralogs/P1.hmm};
results/2026-09-29_aligner_default/.
