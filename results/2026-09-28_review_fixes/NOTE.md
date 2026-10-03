# Review fixes F1, F2, F5 (part), F7, and the A-subloci merge (2026-09-28)

Branch `review-fixes` from PR #9 ad9be39; frozen worktree `run-ca3045c`.
Review: results/2026-09-28_fable_review/README.md.

## Commits
- 4f37cb7  F1: the relaxed pass is gated on strict results that SURVIVE the
  modelled-gene bar and the flank-carried rule (both now run before the gate
  and again on the relaxed calls).
- 558dc1c  F2: an admitted cluster with a modelled gene that falls below the
  fraction floor is listed in `suppressed_loci` (`below_fraction_floor`, with
  fraction, floor, best identity) and in diagnostics, unless a later pass
  reports it. New `detect/record_selfcall.py` and
  `scripts/check_record_selfcall.py` verify each curated record's source genome
  is called at the record coordinates.
- ef6a4bd  F5 (general part): suppressed loci carry `idiomorph_classifier`.
- 18ed551  F7: a merged call is unverified if any member is (reasons
  combined) and takes the best member confidence only when all members carry
  the same idiomorph; otherwise the primary's.
- f487756  A provisional Aalpha/Abeta separation switch (curator 2026-09-28,
  option b) -- superseded by:
- ca3045c  Curator option (a) after the subloci literature review: the
  roster switch is removed (the code keeps `merge_separately`, default off =
  merge); a merged call lists each contributing family under `subloci`
  (label, generic, idiomorph, coordinates, genes, genes missing, completeness).
  Not done (separate later task): grouping distant subloci by conserved flanks
  (mip/beta-fg; S. commune Aalpha-Abeta ~450-550 kb apart).

## F5 spec (for the guard on branch curation-umbelopsis; not implemented here)
`secondary_undetermined.py` should accept as anchor only a same-family call
with `idiomorph_classifier.input == "model"` (not `hsp_fragment`, not a
flank-carried `low` call) and confidence >= the call it would withhold.
Withheld loci now carry the classifier block (ef6a4bd), so the guard's
withheld margins become auditable once curation-umbelopsis is rebased.

## Record self-call check (curation-umbelopsis db, 2e9aa97 reports)
selfcall_curation_umbelopsis_2e9aa97.tsv: called 2 (both Umbelopsis
records), missed 1 (13706_nrrl-2496_MAT_Plus, S. racemosum -- the F2 case;
that run predates the F2 report row), no_report 12, no_assembly 103.
Limitation: 103 of 118 records name no assembly accession anywhere in their
metadata, so the check cannot place them; the schema needs an assembly field.

## Measurements
(see RESULTS below, filled from compare_output.txt)

- f25cf70  Found while measuring F1: relaxed-pass calls ignored the HMM
  classifier (labelled by the old bitscore vote with no classifier block,
  although a verdict is computed for every cluster). They now use it exactly
  as the strict path does.

## RESULTS (compare_output.txt; jobs.txt)

Zygo 23 (ca3045c and f25cf70): scaffold 23/23 locus, 23/23 idiomorph;
contig 23/23, 23/23.

Mucoromycota 293 genomes, f25cf70 vs results/2026-09-27_split_locus/
Mucoromycota_4174440: genomes called 247 -> 258; calls gained 17, lost 0,
changed 13.
- All 17 gains are relaxed-pass calls (medium, 2 modelled genes) in genomes
  whose strict candidates were all withheld (F1). Classifier labels: 9
  undetermined, 5 Plus, 3 Minus. Six genomes get TWO relaxed calls
  (Syncephalastrum spp. x5 incl. PYS2702/2703, S. contaminatum x2; Circinella
  minor; C. umbellata); in each pair one call is undetermined. They look like
  HMG-paralog pairs; core identities are 27-43%. The secondary-undetermined
  guard (curation-umbelopsis) would withhold the undetermined second calls.
- 12 changes are the split-locus calls (11 GL R. arrhizus/delemar + R.
  microsporus 56028): same locus and label, now reported by the relaxed pass
  (sexP + btbA, 2 modelled) at medium instead of split-locus low. The
  split-locus rule no longer fires on them because the relaxed pass reports
  the family first.
- 1 change: Rhizomucor pusillus FCH_5_7 FWWN02000322.1 Minus -> undetermined
  (the relaxed call now takes the classifier verdict).
- F1 cases on this (PR #9) code: A. glauca and R. pusillus FCH_5_7 keep their
  calls (they were lost only on curation-umbelopsis, whose records created the
  blocking noise; re-measure there after rebase).
- suppressed_loci: below_fraction_floor 379 (new, F2), modelled_gene_bar 795,
  flank_carried 12; 1,171 suppressed loci carry a classifier block (F5).

Agaricales CAAX panel (128 genomes), ca3045c + basidio-anchors db vs
results/2026-09-27_locus_merge/agaricales_panel: calls 276 -> 281; merges
unchanged (Aalpha+HD 26, Balpha+Bbeta+PR 5, Bbeta+PR 2); all 33 merged calls
carry `subloci`; unverified calls 29 -> 29; no confidence change on any
matched call (F7 changed no label here). The 5 new calls are relaxed PR calls
(receptor + CAAX precursor) in Pholiotina rufidispora (x2) and Megacollybia
platyphylla (x3), families not on the CAAX-unverified list (review F4). The
panel has no Abeta calls, so A-sublocus separation could not be exercised on
it; the switch is off in the roster.
