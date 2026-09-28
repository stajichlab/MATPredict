# Locus merge and the CAAX unverified label (2026-09-27)

Code: branch locus-merge (a2fe24e merge, 4bf4672 CAAX label) on polish-scope-cuts
44567fc. Measured from a throwaway branch merge-measure = basidio-anchors
(f017ecb, the signed-off PR records) rebased onto locus-merge; frozen worktree
run-merge-measure (88771e5). Baseline: results/2026-09-27_caax_precursor/on_*
(c8412b7 + the same records). Only the two rulings differ.

## Rule 1: same-locus calls reported once
Roster `merge_group` (+ `merge_generic` for the catch-all family):
B = PR (generic), Balpha, Bbeta; A = HD (generic), Aalpha, Abeta. Merge needs
same contig, >= 50% overlap of the shorter span, compatible idiomorphs
(undetermined matches anything). Primary = the specific family; evidence,
genes, span, records unioned; confidence = best member; `merged_from` lists
each member. HD never merges with PR (different groups).

Agaricales panel (125 genomes in both runs): calls 303 -> 268; merges
Aalpha+HD 25, Balpha+Bbeta+PR 4, Bbeta+PR 2. No old locus lost its covering
call (0). Schizophyllum Schco3, H4-8 (GCA_019143615.1) and two other strains
now report one B locus (Balpha+Bbeta+PR); Coprinopsis unchanged apart from its
A locus (Aalpha+HD). Polyporales/Russulales controls (37): 84 -> 84, no merges.
Offline replay: full Basidiomycota run 3,958 -> 3,733 calls (Aalpha+HD 181,
Balpha+Bbeta 44); Ascomycota cap6 611 -> 611.

Ascomycota MAT/MATtub/MATyl/MATsc/MTL duplicates: NOT grouped. cap6 has 80
overlapping cross-family pairs, 63 with differing labels -- mostly the same
idiomorph in different vocabularies (A vs MATa, MAT1-2 vs A). Merging them
needs an idiomorph vocabulary map; not done. Genuine double calls on separate
contigs (e.g. Cercospora kikuchii) never merge.

## Rule 2: CAAX-dependent PR calls in low-enrichment families are unverified
db/caax_unverified_taxa.yml: Agrocybaceae 3710399, Mycenaceae 2024004,
Physalacriaceae 862241, Galerinaceae 3710460 (NCBI family ranks). A PR call
admitted only by counting a caax_scan precursor is labelled; confidence
unchanged. Agaricales panel: 29 calls in 21 genomes labelled (Agrocybaceae 16,
Mycenaceae 6, Physalacriaceae 5, Galerinaceae 2) = every CAAX-gained call in
those families. rollout_summary unverified_calls: 21 genomes, 29 calls.

## Regression
Zygo 23: scaffold 23/23 locus, 23/23 idiomorph; contig 23/23, 23/23.
Tests: 851 pass (840 after rule 1).
