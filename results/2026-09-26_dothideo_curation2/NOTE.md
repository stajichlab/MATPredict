# Dothideomycetes curation, round 2 (2026-09-26)

Branch `curation-mucor-dothideo`, records commit `52cd292`. Before arm
`d18ec5a` (round-1 records signed off; PR #9 code at 90db4af). Same code in
both arms; only the database differs. 100-genome list from
`results/2026-09-26_dothideo_curation/dothideomycetes.tsv`. 4 jobs on `short`,
18-21 min each, 0 failures, 0 suppressed.

## Records added (8, all tier 1 by the database's working bar, pending sign-off)

| record | deposit | genes | source | evidence stated in abstract |
|---|---|---|---|---|
| 5022_unknown-1_MAT_MAT1-1 | AY174048.1 | MAT1-1-1 | Cozijnsen & Howlett 2003, PMID 12679880 | targeted sequencing, transcripts; no cross described |
| 5022_unknown-2_MAT_MAT1-2 | AY174049.1 | MAT1-2-1 | same | same |
| 13684_sn435pl98_MAT_MAT1-1 | AY212018.1 | MAT1-1-1 | Bennett et al. 2003, PMID 12948511 | isolates "known to mate" |
| 13684_sn436ga98_MAT_MAT1-2 | AY212019.1 | MAT1-2-1 | same | same |
| 5499_alenya-b_MAT_MAT1-1 | DQ659350.2 | MAT1-1-1 | Stergiopoulos et al. 2007, PMID 17178244 | targeted sequencing; species presumed asexual |
| 5499_imi-day9-054980_MAT_MAT1-2 | DQ659351.2 | MAT1-2-1 | same | same |
| 1873960_unknown-1_MAT_MAT1-1 | DQ787015.1 | MAT1-1-1 | Conde-Ferraez et al. 2007, PMID 20507483 | targeted sequencing; no cross described |
| 1873960_unknown-2_MAT_MAT1-2 | DQ787016.1 | APN2, MAT1-2-1 | same | same |

All 9 proteins validate at 100% identity and coverage. Partial L. maculans
DNA lyase (codon_start=2, 3'-partial) left out. Annotation errors C8-C10.

## Before/after

| | before d18ec5a | after 52cd292 |
|---|---|---|
| genomes called | 89 | 87 |
| loci high / medium | 47 / 43 | 44 / 46 |
| genomes with both idiomorphs called | 1 | 3 |

8 genomes changed (`compare_output.txt`):

- **2 corrections.** P. nodorum SnOre 11-1 and Mur_S3: MAT1-1 (margin 1.07)
  -> MAT1-2 (margin 343.9). Direct tblastn: P. nodorum MAT1-2-1 hits both
  genomes at 98.9%; P. nodorum MAT1-1-1 has no hit at E <= 1e-20. The old
  MAT1-1 calls were wrong.
- **2 losses caused by the polish cap, not by the records.** Z. brevis
  (target cluster 95.6% identity) and Acidiella bohemica (81.8%). The new
  references admit more clusters (32 -> 40; 58 -> 65), and the genes-first cap
  of 6 now skips the true locus, which is admitted but `polish_capped`.
- **2 new double calls, unverified.** Cercospora kikuchii gains a MAT1-2
  medium call (contig VTAY01000089.1) next to its MAT1-1 call; Nothopassalora
  personata gains a 940 bp MAT1-2 partial_locus next to its MAT1-1 call.
  Both species are expected heterothallic; the second calls are probably
  spurious.
- **1 flip at a tiny margin.** Aureobasidium pullulans: MAT1-2 high (margin
  5.8) -> MAT1-1 medium (2.3). Both MAT1-1-1 and MAT1-2-1 are in the locus;
  Aureobasidium is homothallic (Gostincar et al. 2014), so neither label fits.
- **1 confidence drop.** Neodevriesia sp. NU419: MAT1-2 high -> medium.

## The 9 genomes uncalled in the round-1 pilot: none fixed

All fail the modelled-gene bar (0-1 genes modelled).

| genome | cause |
|---|---|
| Tubeufia hainanensis GCA_975972555.1 | not a genome: 0.1 Mb MAG bin, 4 contigs (suppress candidate) |
| Cladosporiales sp. GCA_024708745.1 | assembly: 12,770 contigs, N50 3.9 kb |
| Aulographum hederae GCA_010015705.1 | its best cluster (83.6%) is polish-capped; also 0.43 Mb of N |
| Bathelium albidoporum GCA_021031095.1 | reference gap (best cluster 43% identity, capped) |
| Microcyclospora pomicola, Piedraia hortae, Vermiconidia calcicola | best cluster polished (64-80% identity) but only 1 gene models: reference divergence (Capnodiales, Extremaceae uncurated) |
| Lineolata rhizophorae, Lojkania enalia | best cluster polished (82-84%) but 0 genes model; 0.59 Mb and 3.27 Mb of N |

Two genomes called in round 1 (`dca6ccc`) were already lost on the newer PR #9
code before these records: Lopnu1 (withheld by the flank-carried rule;
suppressed_flank_carried = 1) and Bauco1 (best cluster polish-capped).

## Limits

- Evidence for the 8 records is from abstracts; two papers were not readable.
- "Correction" is shown only for the 2 P. nodorum genomes. The double calls
  and the Aureobasidium flip were not checked against any truth.
