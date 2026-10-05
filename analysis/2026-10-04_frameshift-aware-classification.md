# Frameshift-aware classification of exonerate models

Status: awaiting curator review (code change).

## Question
A Nanopore homopolymer frameshift truncated the Mucor griseocyanus CBS 116.08
sexM model, and the call fell to undetermined (notable finding 026). Can the
classifier score the frame-corrected protein instead, and what else changes?

## Ruling
J. Stajich, 2026-10-04: implement the proposed fix.

## Method
- Exonerate reports frameshifts as an exon attribute (`frameshifts N`) and its
  aligned target blocks on the `similarity` line (`Align t q len`). Block
  coordinates measured on exonerate 2.4.0 with the M. griseocyanus sexM on both
  strands: plus [t, t+len-1], minus [t-len, t-1].
- `PolishModel` gains `frameshifts` and `cds_blocks`; `_translate_model` uses
  the joined blocks (frameshift bases and target-only insertions left out) when
  the model has frameshifts and the result has no internal stop.
- Reports: gene evidence `frameshifts: N`; classifier verdict
  `frameshift_corrected: true` (only when present).
- Code: branch `frameshift-aware` (84b6b55, on main after #25); tests 1022
  passed (3 new).

## Results
- M. griseocyanus CBS 116.08: undetermined -> **Minus (high)**, scores
  Minus 154.2 vs Plus 43.5 (margin 110.7); sexM `frameshifts: 5`.
- 288 BFD Mucoromycotina (campaign mode) vs the same code without the fix
  (356b9f3): 8 genomes have a frameshift-corrected call; 1 call changes (the
  one above); the other 7 keep their idiomorph
  (`results/2026-10-04_mucoro_bfd_frameshift/compare.txt`).
- Regression panel + Zygo 23 vs v0.6.1: identical to the #25 diff (no
  additional change); Zygo 23/23 on both inputs
  (`results/2026-10-04_regression_frameshift/`).
- Ascomycota and Basidiomycota: 0 changes (their idiomorph comes from the
  search path, not the HMM classifier).

## Limits
- A frameshift is corrected, not proven to be an assembly error; the flag says
  "possible assembly error".
- miniprot models are not frameshift-corrected (miniprot reports frameshifts
  differently; not needed for this case).

## Decision
Pending curator review.
