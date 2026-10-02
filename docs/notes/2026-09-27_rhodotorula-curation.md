# Rhodotorula P/R and HD MAT records (2026-09-27)

Curator request: curate the Rhodotorula MAT loci, including HD, from the
group's manuscript (Liu, Tsai, Coelho, ... Stajich; bioRxiv
10.1101/2025.09.11.675505, submitted, not yet peer-reviewed; data on Zenodo
10.5281/zenodo.18230104). The manuscript itself is not in the repository.

## Records (tier 2, genome-derived, pending curator sign-off)

All genes are recorded by coordinates on public contigs; none has a protein
accession in the public assemblies. Each is either the group's own gene model
or, where that model did not translate on the public contig, a miniprot 0.18
re-model (named in the record). All 21 translations were checked: no internal
stop, and `curate-db build-gff` derives the same proteins.

| record | assembly | genes |
|---|---|---|
| 5286_cbs-14_redPR_A1 | GCA_921037615.3 | STE3a1 |
| 5286_nbrc-0880_redPR_A2 | GCA_000988875.2 | STE3a2, STE20 |
| 5535_cbs-20_redPR_A1 | GCA_920103745.3 | STE3a1, STE20 (miniprot) |
| 5537_jy1105_redPR_A1 | GCA_024748845.1 | STE3a1, STE20 |
| 86839_nrrl-y-7192_redPR_A2 | GCA_060419215.1 | STE3a2, STE20 |
| 29898_jj10-1_redPR_A2 | GCA_026119225.1 | STE3a2, STE20 |
| 5286_cbs-14_redHD_B41 | GCA_921037615.3 | HD1, HD2 (miniprot) |
| 5286_nbrc-0880_redHD_B39 | GCA_000988875.2 | HD1, HD2 |
| 5535_cbs-20_redHD_B63 | GCA_920103745.3 | HD1, HD2 |
| 5537_jy1105_redHD_B8 | GCA_024748845.1 | HD1, HD2 (miniprot) |
| 86836_ls11_redHD_B74 | GCA_002917965.1 | HD1, HD2 |

- Both P/R alleles in each of the manuscript's three clades (A: JY1105 A1,
  Y-7192 A2; B: CBS 14 A1, NBRC 0880 A2; C: CBS 20 A1, JJ10.1 A2).
- Every new STE3 hits the curated Coelho 2011 receptor of the same allele
  first (55-78% identity). HD1 hits the rust bW (HD1) and HD2 the rust bE
  (HD2) first, which fits the HD1/HD2 classes.
- JJ10.1 is *Rhodotorula graminis* in NCBI; the manuscript places it in
  Rhodotorula sp. (clade I, R. aff. babjevae). The record uses the NCBI taxid
  and says so.
- Left out (B5.6): LS11 STE3.A2 and Y-7192 HD2 (frameshifts in the public
  assemblies), CBS 14 STE20 (internal stop). R. mucilaginosa Y-2510 and the
  other newly sequenced reference strains have no public assembly yet
  (only R. sphaerocarpa Y-7192 was released, 2026-08-25).

## Roster changes (db/Basidiomycota/order.yml)

- `redPR`: STE20 declared as an optional flank, `idiomorph_informative:
  false`, and `exclude_from_search: true`, following the Tremellales precedent
  (2026-09-21: STE20 is a large ubiquitous kinase family; 187 of 272
  genome-wide HSPs on Cryptococcus JEC21). Cluster gap unchanged (20 kb).
- New `redHD` family: Sporidiobolales, pattern `^B[0-9]+$` (the manuscript's
  HD allele names), HD1 + HD2 core, 10 kb gap. Fills the queued
  Sporidiobolales HD gap.
- Pheromone precursors (RHA) are not curated; they belong to the queued
  receptor/pheromone work.

## Validation (results/2026-09-27_rhodotorula_curation/)

221 genomes: all 216 BFD Rhodotorula (suppressed skipped) plus the 15
Sporidiobolales of the 2026-09-26 test. Before = `ac33880`, after =
`02434ee`, only the database differs.

| | before | after |
|---|---|---|
| P/R called | 215 | 216 |
| HD called | 0 | 217 |
| median wall per genome | 18 s | 106 s |

- **HD and P/R never land on the same contig** (0 of the genomes with both),
  as the manuscript reports for the genus (physically separate loci).
- **Known answers:** 67 BFD genomes have a strain allele assignment in the
  manuscript's tables. P/R agrees for 66/67 in both arms (held out, i.e.
  excluding the 7 record strains: 61/62). The one miss is a manuscript
  A1/A2 hybrid (R. mucilaginosa RIT389) that gets no P/R call; it gets one
  HD call. The other hybrid in the set (R. toruloides CCT 0783) is called
  A1+A2 with two HD calls, matching the manuscript. HD is called in 67/67
  (HD allele identity is not callable by homology).
- **Clade A A2 dominance:** R. mucilaginosa 117 A2, 6 A1, 2 uncalled.
- HD uncalled in 4 genomes: Sporobolomyces pararoseus (outside Rhodotorula),
  one R. mucilaginosa, one R. toruloides and an R. graminis MAG.
- P/R call rate was already 215/221 with the Coelho 2010/2011 records; the
  gain from the new P/R records is 1 genome. The gain from this work is HD.
- **Cost:** adding the HD search raises the median wall 5.9x (18 -> 106 s;
  range seen 60-410 s). Most redHD clusters are genome-wide homeodomain hits
  that are not admitted (e.g. 71 of 75 in LS11).
