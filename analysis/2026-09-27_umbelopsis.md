# Umbelopsis: lost calls, tier-2 records, the guard, and Circinella
Status: decided (records signed off 2026-09-29; guard dropped); running (P1 update on curation-umbelopsis)

## Question
Why did Umbelopsidaceae drop from 13/14 to 8/14 called, do tier-2 records
recover them, and what happened to Circinella minor?

## Data and code version
- Branch `curation-umbelopsis` (pushed), rebased onto PR #9 076afe4 and then
  9a458de; records `41833_gzumbrama1_MAT_Plus`, `44442_wa0000051536_MAT_Minus`,
  `13706_nrrl-2496_MAT_Plus`.
- Scan: 293 Mucoromycota genomes; Zygo 23 on both inputs.

## Results
- Diagnosis (`results/2026-09-27_umbelopsis_diagnosis/NOTE.md`): the flank
  rule's genome-size-dependent E <= 1e-5 floor withheld 5 loci with the Mucorales
  gene order. Fixed by the 39-bit floor (flank-carried report).
- Guard (`results/2026-09-28_umbelopsis_rebased/NOTE.md`): under the MAT-gene
  gate the secondary-undetermined guard withheld nothing (guard on vs off:
  253 vs 253 genomes, 0 calls changed; Zygo 23/23 both). Dropped (revert).
- Effect of the records: Umbelopsidaceae 13/14 called (was 12, of which 4
  undetermined); all model-typed; 7 Plus, 6 Minus. U. vinacea x2 and WA50703 x2
  now Minus/high; U. nana newly called; AD052 still not called.
- Record self-call: all three records call their own source genome
  (gzUmbRama1 Plus/high; WA0000051536 Minus/high; S. racemosum NRRL 2496
  Plus/medium).
- Circinella minor (`results/2026-09-29_circinella_trace/NOTE.md`): lost
  because the classifier rebuild lowered its sexP score from 112.5 to 85.9
  bits, below the gate, with one flank (rnhA). The locus and model did not
  change. (This supersedes the rebase note's "reshaped cluster" guess.) Its
  tptA/algA/glrA sit on other contigs, the S. racemosum pattern; withholding
  is consistent with the Lichtheimiaceae ruling (a).
- The gate threshold per build (98.9 bits on this branch) changed no call
  (`results/2026-09-29_gate_threshold/NOTE.md`).

## What changed in detection
Records and classifier rebuild on curation-umbelopsis; guard commit reverted.
P1 class and gate threshold being added via `--paralogs-only` / `--gate-only`
(`results/2026-09-29_umbelopsis_p1/`, running at the time of writing).

## Limits
- Four Umbelopsis Minus margins are 31.9-37.2 bits, near the 25-bit floor.
- The Minus record's sexM is unannotated and may be truncated.
- 52 calls shifted span slightly after the records were added; not inspected
  one by one.
- This branch's HMMs were built before builds were deterministic.

## Curator decisions
- Made 2026-09-29: drop the guard; sign off all three records; trace
  Circinella (done).
- Open: merging curation-umbelopsis into PR #9 after the P1 update.

## Files
`results/2026-09-27_umbelopsis_diagnosis/`, `results/2026-09-27_umbelopsis_curation/`,
`results/2026-09-27_mucoro_curation_guard/`, `results/2026-09-28_umbelopsis_rebased/`,
`results/2026-09-29_circinella_trace/`.
