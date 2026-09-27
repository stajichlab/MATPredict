# 002. Rhodotorula P/R and HD loci recovered genome-wide; agreement with an independent study

- **Category:** curation-first
- **Status:** verified (against the study's published-preprint strain assignments)
- **Lineage:** Basidiomycota, Pucciniomycotina, Sporidiobolales

## Summary
Eleven tier-2 records (six P/R, five HD) from public Rhodotorula assemblies let
detection call the HD locus in 217 genomes (0 before). The P/R allele agrees
with the study's independent assignments for 61 of 62 held-out strains, and
HD and P/R never share a contig, consistent with physically separate loci.

## Evidence
- Records (family `redPR`): `5286_cbs-14_redPR_A1` (GCA_921037615.3),
  `5286_nbrc-0880_redPR_A2` (GCA_000988875.2), `5535_cbs-20_redPR_A1`
  (GCA_920103745.3), `5537_jy1105_redPR_A1` (GCA_024748845.1),
  `86839_nrrl-y-7192_redPR_A2` (GCA_060419215.1), `29898_jj10-1_redPR_A2`
  (GCA_026119225.1).
- Records (new family `redHD`, HD1+HD2): `5286_cbs-14_redHD_B41`,
  `5286_nbrc-0880_redHD_B39`, `5535_cbs-20_redHD_B63`, `5537_jy1105_redHD_B8`,
  `86836_ls11_redHD_B74` (GCA_002917965.1).
- New STE3 proteins hit the same-allele Coelho et al. 2011 receptor first
  (55-78% identity); HD1 hits rust bW (HD1) and HD2 rust bE (HD2) first.
- Validation on 221 genomes (216 BFD Rhodotorula + 15 Sporidiobolales), only the
  database differing: HD 0 -> 217 called; P/R 215 -> 216; HD and P/R on the same
  contig in 0 genomes. Source: `docs/notes/2026-09-27_rhodotorula-curation.md`
  (curation-puccinio); tables in `results/2026-09-27_rhodotorula_curation/`
  (main checkout; answer table kept out of the repository).
- Known answers: P/R agrees 66/67 (61/62 excluding the 7 record strains). The
  miss is a reported A1/A2 hybrid (entry 010).
- *R. mucilaginosa*: 117 A2, 6 A1, 2 uncalled (clade A dominated by A2, as the
  preprint reports).

## Method that found it
Coordinates-on-public-contig records (three genes re-modelled with miniprot
0.18 where the group's model did not translate); detection before/after with
frozen worktrees `run-ac33880` and `run-02434ee` (curation-puccinio `02434ee`).
Source study: Liu et al., bioRxiv doi:10.1101/2025.09.11.675505; data Zenodo
doi:10.5281/zenodo.18230104. Curator sign-off 2026-09-27 (`16cfcbe`).

## Verification done / still open
- Done: agreement with the study's strain assignments (61/62 held out).
- Open: HD allele identity is not callable by homology (HD called 67/67 where
  answers exist, but not typed); the newly sequenced reference strains
  (e.g. *R. mucilaginosa* Y-2510) have no public assembly yet.

## Limits
- The comparison set is a preprint, not yet peer-reviewed.
- The HD search costs 5.9x runtime (median 18 -> 106 s per genome); most of it
  is the genome-wide tblastn (`results/2026-09-27_hd_prescreen/NOTE.md`).
- RHA pheromone genes are not curated (entry 003).

## Related
Entry 010 (hybrids); annotation report B5.6 (frameshifts in public assemblies).
