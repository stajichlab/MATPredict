# curation-umbelopsis rebased onto PR #9 076afe4 (2026-09-28)

## What was done

- Rebased `curation-umbelopsis` onto 076afe4 (the review fixes, the MAT-gene
  gate, F6, the all-CAAX unverified label and the assembly field). A backup of
  the pre-rebase branch is `backup/curation-umbelopsis-pre-076afe4`.
- One conflict in `pipeline.py`: the MAT-gene gate import and the
  below-fraction-floor block (from PR #9) against the secondary-undetermined
  guard (this branch). Both were kept. The guard now runs before the floor
  block, so the floor block's "reported after all" check sees final results.
- Classifier rebuilt with `scripts/build_idiomorph_hmms.py`: LOO 88/88 correct,
  worst correct margin 13.7 bits (was 13.5). sexP n=78 (27 genera), sexM n=10
  (8 genera). Tests: 906 passed.

## Is the secondary-undetermined guard still needed? No.

Two frozen worktrees, same code except the guard commit reverted:
`run-667252a` (guard on) and `run-umb-noguard` (guard reverted, 9aa90a2).

| | guard on | guard off |
|---|---|---|
| Mucoromycota genomes called (293) | 253 | 253 |
| calls gained / lost / changed between arms | 0 / 0 / 0 | |
| genomes with an undetermined call beside a determined one | 0 | 0 |
| Zygo 23, scaffolds | 23/23 locus, 23/23 idiomorph | 23/23, 23/23 |
| Zygo 23, contigs | 23/23, 23/23 | 23/23, 23/23 |

The MAT-gene gate (>= 100 bits, or >= 2 flanks modelled at >= 40%) already
removes every spurious second call the guard targeted. The guard withholds
nothing. Proposal: drop it (revert f194206 on this branch; the reverted state is
measured as `umb-noguard`).

## Effect of the records (guard arm vs the 076afe4 baseline scan)

253 genomes called in both. 52 calls keep the same locus but a shifted span
(the new flank references widen clusters); same label and confidence unless
listed below. 1 new locus, 1 lost.

Umbelopsidaceae (14 genomes):

| Genome | Before | After |
|---|---|---|
| U. isabellina B7317, MPG-14A | Plus/low (fragment) | Plus/high (model) |
| U. isabellina WA0000067209 | Plus/low (fragment) | Plus/medium (model) |
| U. isabellina NBRC 7884 | Plus/high | Plus/high |
| U. sp. M5902 | Plus/low (fragment) | Plus/high (model) |
| U. vinacea WA0000051536 (record source) | undetermined/low | Minus/high (margin 99.5) |
| U. vinacea gzUmbVina2 | undetermined/low | Minus/high (36.5) |
| U. sp. WA50703 x2 | undetermined/low | Minus/high (37.2) |
| U. nana NBRC 117090 | not called | Minus/high (31.9), new |
| U. ramanniana gzUmbRama1 (record source) | Plus/high | Plus/high (135.3) |
| U. ramanniana AG (Umbra1) | Plus/high | Plus/high (125.7) |
| U. sp. PMI_123 | Plus/high | Plus/high |
| U. sp. AD052 | not called | not called |

13 of 14 called (was 12, of which 4 were undetermined). All calls are now
model-typed; 7 Plus, 6 Minus (5 of the Minus margins are 31.9-99.5).

Syncephalastraceae: every call keeps its label and confidence; margins rise
(e.g. S. racemosum B6101 125.4 -> 186.3; S. monosporum B8922 Minus 29.7 -> 34.5).

## Record self-call (scripts/check_record_selfcall.py)

| Record | Source genome | Verdict |
|---|---|---|
| 41833_gzumbrama1_MAT_Plus | GCA_977110945.1 | called, Plus/high |
| 44442_wa0000051536_MAT_Minus | GCA_016758895.1 | called, Minus/high |
| 13706_nrrl-2496_MAT_Plus | GCA_002105135.1 | called, Plus/medium |

## Previously flagged side effects

- R. pusillus FCH_5_7: not called in either the baseline or this branch
  (Lichtheimiaceae; the gate needs flanks they lack — out of scope, research
  follow-up).
- A. glauca: Minus/medium in both (margin 106.1 -> 93.3).
- M. griseocyanus CBS 116.08: Minus/high in both (margin 28.5 -> 26.3;
  still just above 25).
- gzUmbRama1 and Umbra1: high (F6 fixed the distant rnhA cap).
- R. stolonifer PRFJ01/PRFJ02: Minus/high in both baseline and branch.
- New loss: Circinella minor GCA_016758965.1 (Plus/medium in baseline) is
  withheld by the MAT-gene gate (margin 53.6, below 100 bits, insufficient
  flanks) once the new flank references reshape its cluster. Not diagnosed.

## Limits

- Umbelopsis Minus margins for 4 genomes are 31.9-37.2, close to the 25-bit
  typing floor; the Minus record's sexM is unannotated and may be truncated.
- Span shifts in 52 calls were not inspected one by one.
