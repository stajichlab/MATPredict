# curation-umbelopsis rebased onto PR #9 R4, with the P1 paralog class

## Question
Does curation-umbelopsis (the Umbelopsis Plus/Minus and S. racemosum NRRL 2496
records, their classifier) still work on top of PR #9 da35821 (R4 P1 paralog
class, per-build gate threshold, two_idiomorphs), and what does it change
against PR #9 alone?

## Data and code version
- Branch curation-umbelopsis at 92eb2a6 (rebased onto da35821); frozen worktree
  .claude/worktrees/run-92eb2a6.
- Baseline: PR #9 R4 runs, code 41bd471 (results/2026-09-29_r4_paralog/).
- Inputs: Mucoromycota 293 (results/2026-09-26_early_diverging/lists/Mucoromycota.tsv,
  1 suppressed genome skipped), Zygo 23 (both inputs), LCG Mucoromycotina 621
  (results/2026-09-28_lcg_holdout/lists2), Jena 64 scaffolds.
- SLURM jobs 29213213-29213217 (jobs.txt).

## Method
- Rebase: the guard commit and its revert were skipped (same end state; the
  guard file is absent). Classifier-file conflicts kept this branch's MAT HMMs:
  sexM.hmm, sexP.hmm and training_extra.faa are byte-identical to 6e5034b
  (sha256 9124f62b4b0b, 87660ee54479, fca74103a410). No full rebuild.
- P1 added with `build_idiomorph_hmms.py --paralogs-only`: 0/88 MAT training
  proteins classed as P1; self 1/1. The P1 HMM equals PR #9's except the header
  date and checksum lines.
- Gate refreshed with `--gate-only`: min_score 98.9 bits; 9/189 paralog
  negatives and 76/88 held-out MAT proteins at or above.
- compare.py compares each arm with the R4 run; changes.tsv, compare_output.txt,
  umbelopsidaceae.tsv, selfcall.tsv.

## Results
- Tests: 930 passed.
- Zygo 23: 23/23 locus and 23/23 idiomorph on scaffolds and on contigs.
- Record self-check (scripts/check_record_selfcall.py, selfcall.tsv):
  41833 gzUmbRama1 Plus/high; 44442 WA0000051536 Minus/high; 13706 S. racemosum
  NRRL 2496 Plus/medium. All three call their own source genome.
- Mucoromycota 293: called 253 -> 253; 1 gained, 1 lost, 8 label/confidence
  changes, 44 span-only changes (same place, same label and confidence).
  - Gained: U. nana Minus/high.
  - Lost: Circinella minor GCA_016758965.1 (gate; only call).
  - Changed: 4 Umbelopsis undetermined/low -> Minus/high (U. vinacea x2,
    WA50703 x2); 4 Umbelopsis Plus/low -> Plus/high or medium.
- Umbelopsidaceae: 12/14 -> 13/14 called; all 13 are model calls; AD052 still
  uncalled. Margins rise (e.g. gzUmbRama1 56.1 -> 135.3). Thin Minus margins:
  WA50703 x2 37.2, gzUmbVina2 36.5, U. nana 31.9.
- LCG 621: called 536 -> 535; 3 gained, 5 lost, 1 changed, 102 span-only.
  - Gained: Mucor ramannianus NRRL A-21216 Minus/high (an Umbelopsis by modern
    name), Mucor sp. NRRL 1454 Minus/high, S. racemosum NRRL 2495 Minus/medium.
  - Lost (all Plus/medium, now withheld by the MAT-gene gate + fraction floor
    unless noted): Circinella minor CBS 143.56 and NRRL 1365, C. umbellata
    NRRL 2417, Rhizopus microsporus NRRL A-17693, and Mucor pusillus NRRL A-13674
    (modelled_gene_bar). 4 of the 5 were the genome's only call.
  - Cause of the gate losses: the branch classifier scores the same sexP+rnhA
    protein 85.3 bits (R4: 112.0), below the 98.9 gate, with one flank. R.
    microsporus A-17693 has the identical protein score profile as C. minor
    NRRL 1365 (112.0/30.4 -> 85.3/32.5); its species label may be wrong (not
    checked).
  - Changed: Umbelopsis ovata NRRL 13127T undetermined/low -> Minus/high.
- Jena 64: called 61 -> 61; 1 change (CBS206_69 undetermined/low -> Minus/low).
- P1 withholdings: Mucoromycota 13 -> 23, LCG 33 -> 39, Jena 8 -> 8. Only one
  withheld P1 locus leaves a genome uncalled (LCG Mucor laxorrhizus NRRL 2814;
  also uncalled in the R4 run? not checked).
- Clean LCG strain labels (n=17; disputed, misidentified, leaked excluded):
  10 agree / 3 disagree / 4 uncalled, unchanged from R4.

## What changed in detection
Commit 92eb2a6 (manifest: P1 class and 98.9-bit gate on the branch HMMs). No
code change on the branch beyond PR #9.

## Limits
- The branch HMMs were built before R4 and are not bit-reproducible on rebuild.
- The Circinella/Rhizomucor losses fall under the Lichtheimiaceae known gap
  (ruling a); R. microsporus A-17693 does not, unless mislabelled.
- M. pusillus NRRL A-13674 lost a strong call (Plus 199.7 bits in R4) to the
  modelled-gene bar; not diagnosed.
- Span-only changes (146 in total) were not checked one by one.

## Curator decisions
- Made: records signed off; guard dropped; P1 via --paralogs-only; gate via --gate-only.
- Open: whether to merge curation-umbelopsis into PR #9 given the LCG trade
  (+3 calls, 5 Umbelopsis-group calls upgraded, 5 Circinella-type/Rhizomucor
  calls lost); diagnose M. pusillus A-13674 and R. microsporus A-17693.

## Files
compare.py, compare_output.txt, changes.tsv, umbelopsidaceae.tsv, selfcall.tsv,
jobs.txt, run_jena.slurm, run_lcg_chunk.slurm; runs in Mucoromycota_92eb2a6/,
zygo23_92eb2a6/, lcg_runs/, jena/runs_scaffolds/.
