# Umbelopsis: lost calls, tier-2 records and the guard
Status: running (re-measure on 076afe4 pending; records await sign-off)

## Question
Why did Umbelopsidaceae drop from 13/14 to 8/14 called, and do tier-2
records recover them?

## Results
- Diagnosis (`results/2026-09-27_umbelopsis_diagnosis/NOTE.md`): the flank
  rule's genome-size-dependent E<=1e-5 floor withheld 5 loci with the Mucorales
  gene order (core HMG fragments 33–44 bits).
- Fixed by the 39-bit floor (see flank-carried report).
- Records (branch curation-umbelopsis, not pushed): Plus
  `41833_gzumbrama1_MAT_Plus`, Minus `44442_wa0000051536_MAT_Minus`.
  With them all 14/14 Umbelopsidaceae were called, plus 13 spurious second
  calls (`results/2026-09-27_umbelopsis_curation/NOTE.md`).
- Guard (bbcb524) withheld 17 spurious second calls
  (`results/2026-09-27_mucoro_curation_guard/NOTE.md`); S. racemosum record
  `13706_nrrl-2496_MAT_Plus` added; Lichtheimiaceae failed the gene-order test.
- Re-measure after rebase onto 076afe4: running
  (`results/2026-09-28_umbelopsis_rebased/`).

## Curator decisions
Open: sign-off of the three records after the re-measure; whether the guard
is still needed under the MAT-gene gate.
