# 006. Cercospora kikuchii assembly carries MAT1-1-1 and MAT1-2-1 at distinct loci

- **Category:** homothallism candidate
- **Status:** candidate
- **Lineage:** Ascomycota, Dothideomycetes, Mycosphaerellales

## Summary
One C. kikuchii assembly gives two MAT calls on two contigs, and both
idiomorphs' core genes are modelled. The second locus carries MAT1-2-1 and a
MAT1-1-1 copy. The cause is not resolved.

## Evidence
- GCA_009193115.1. `VTAY01000105.1`: MAT1-1 call, MAT1-1-1 64.0% (modelled),
  COX13 72.6%. `VTAY01000089.1`: MAT1-2 call, MAT1-2-1 70.8% (modelled),
  MAT1-1-1 51.2% (modelled), MAT1-1-3 34.4% and APN2 29.0% (unmodelled).
- Source: `results/2026-09-27_dothideo_recheck/NOTE.md` (committed on
  curation-mucor-dothideo `66796e6`); first seen in
  `results/2026-09-26_dothideo_curation2/`.

## Method that found it
Detection with the curated Dothideomycetes records (Cochliobolus, Zymoseptoria,
Parastagonospora, Leptosphaeria, Pseudocercospora) on a 100-genome pilot.

## Verification done / still open
- Done: persisted after the final code (tier rule, cap rank); the parallel
  double call in Nothopassalora personata was an artefact over COX13 and
  disappeared.
- Open: homothallic arrangement vs mixed assembly vs paralog. Check with reads
  and literature.

## Limits
One assembly; no reads examined.

## Related
Entry 005; curation note `docs/notes/2026-09-26_mucorales-dothideomycetes-curation.md` (curation-mucor-dothideo).
