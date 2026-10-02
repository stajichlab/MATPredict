# Mucoromycota-group scan, chytrid control and flank-synteny check
Status: decided

## Question
How does detection perform across Mucoromycota, Mortierellomycota,
Kickxellomycota and a chytrid negative control; are the non-Mucorales calls
real MAT loci?

## Data and code version
- Code: 634dda4 (frozen worktree run-634dda4); `--phylum Mucoromycota` for the
  non-Mucoromycota phyla; chytrids with default routing.
- Inputs: BFD genomes; lists in `results/2026-09-26_early_diverging/lists/`.

## Method
1. Detection panel per lineage.
2. Flank-ortholog synteny check: confirm tptA/rnhA orthologs by best reciprocal
   match to Mucorales proteomes; test whether an HMG hit lies between them.

## Results
Scan (`results/2026-09-26_early_diverging/summary.txt`):

| Lineage | Called |
|---|---|
| Mucoromycota | 227/293 |
| Mortierellomycota | 6/100 |
| Kickxellomycota | 7/190 |
| Chytrid control | 0/19 (6/25 timed out under exhaustive routing) |

Synteny (`results/2026-09-26_flank_synteny_ED/NOTE.md`): tptA and rnhA within
100 kb with an HMG hit between them in 117/293 Mucoromycota, 0/100
Mortierellomycota, 0/190 Kickxellomycota, 0/19 chytrids. None of the 13
Mortierellomycota/Kickxellomycota calls is supported. Lichtheimiaceae (0/30),
Umbelopsidaceae (0/14) and Syncephalastraceae (0/9) also lack the pair.

Superseded: the note's Rhizopus Plus skew was first read as sampling bias. It
is a labelling artefact (btbA vote); see
`2026-09-26_idiomorph-labelling-and-classifier.md`.

## What changed in detection
- `f737863`: calls from out-of-phylum forced routing are `verification: unverified`.
- Not-searched routing replaced exhaustive routing (`3684c60`).

## Limits
A missing tptA–rnhA pair does not prove absence (moved loci are invisible to
this test). Every genome carries 45–69 sexM/sexP-like HMG hits.

## Curator decisions
Made: Mortierellomycota, Kickxellomycota, Zoopagomycota, Glomeromycota,
Blastocladiomycota, Chytridiomycota and Endogonales are discovery-only.

## Files
`results/2026-09-26_early_diverging/`, `results/2026-09-26_flank_synteny_ED/NOTE.md`,
`docs/notes/2026-09-26_early-diverging-scan-and-chytrid-control.md`.
