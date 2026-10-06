# Independent review (Fable model), 2026-09-28

Read-only review of PR #9 (ad9be39) rule changes from 2026-09-26/27 and the
curation-umbelopsis branch (4a58b88). The reviewer ran both test suites
(851 and 836 passed). The full report was delivered in-session; findings are
summarised here.

Verified by the coordinator: F1 (relaxed-pass gate at pipeline.py:2690 runs
before the modelled-gene bar and the flank-carried rule at 2758).

## Findings (severity; status)
- F1 MAJOR, verified: relaxed pass gated before withhold rules; new flank
  references create strict noise clusters that are later withheld but block
  the relaxed pass. Causes the R. pusillus FCH_5_7 and A. glauca losses.
- F2 MAJOR, verified: S. racemosum NRRL 2496 record locus dropped at the 0.5
  fraction floor (sexP+rnhA = 2/5) with no report row. Needs a report row for
  admitted+polished clusters dropped at the floor, and a test that each
  record's own genome is called.
- F3 MAJOR: min_margin 25 is circular (7 cases, weak FastTree) and does not
  separate HMG paralogs on fragment models; the guard patches over it.
- F4 MAJOR: CAAX admission has an unbudgeted false-positive rate (~45 of 118
  gained calls likely non-mating; only 29 labelled). Needs a negative control.
- F5 MAJOR: guard anchor too permissive (any determined call, incl.
  flank-carried or fragment verdicts); suppressed loci omit classifier block.
- F6 MAJOR/MINOR: distant unpolished flank HSPs inside the 50 kb gap cap
  loci at medium (gzUmbRama1, Umbra1). General hazard.
- F7 MINOR: merge keeps only the primary's verification (drops CAAX
  unverified); primary at low takes max confidence. Group A needs sign-off.
- F8 MINOR: split-locus sound; tuned on one Rhizopus batch.
- F9 MINOR: two tests assert roster state, not intent.
- Zygo 23 is saturated (23/23 under every rule) and no longer discriminates.
- Q5: Ascomycota idiomorph vocabularies need a per-locus synonym map.
