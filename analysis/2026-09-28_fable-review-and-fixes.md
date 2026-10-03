# Independent (Fable) review and fixes
Status: decided (fixes landed); some items open

## Results
- Review (`results/2026-09-28_fable_review/README.md`): 6 major findings.
- Fixes (`results/2026-09-28_review_fixes/NOTE.md`, 3aec88b): F1 relaxed gate,
  F2 floor-dropped clusters reported, F5 classifier block on suppressed loci,
  F7 merged verification, subloci evidence, relaxed-pass classifier.
  Mucoromycota 247->258 genomes.
- Zygo 23 saturated (23/23 under every rule); new held-out sets needed.

## Details (added 2026-09-29)
- Review findings (`results/2026-09-28_fable_review/README.md`): F1 relaxed
  pass gated before the withhold rules (verified at pipeline.py); F2 a record's
  own locus dropped silently at the fraction floor; F3 the 25-bit margin is
  circular; F4 CAAX false positives unbudgeted; F5 guard anchor too
  permissive; F6 distant unpolished flank hits cap confidence; F7-F9 minor.
- Review fixes (`results/2026-09-28_review_fixes/NOTE.md`, fast-forward to
  3aec88b; 882 tests): Mucoromycota genomes called 247 -> 258 (17 gained,
  0 lost); 379 floor-dropped clusters now reported; 1,171 suppressed loci carry
  a classifier block; relaxed-pass calls now use the classifier; group A merge
  kept with subloci as structured evidence (see subloci report).
- Follow-up fixes (`results/2026-09-28_next_fixes/NOTE.md`, 076afe4; 900
  tests): see the MAT-gene gate report.

## Curator decisions
Made 2026-09-28: fix order F1, F2, F5, F7, F6, then validation F3/F4, then
Q5/F9. Open: Q5 (Ascomycota idiomorph synonym map), F9 (tests asserting
roster state).
