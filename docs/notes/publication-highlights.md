# Candidate findings to call out in a MATPredict publication

A running list, kept by curator request. Each entry names the evidence and its
limits. Nothing here is a claim until the evidence is re-checked for the paper.

## Wallemia MAT locus detected genome-wide (2026-09-27)

- A putative MAT locus (SXI1, HMG, STE3 with BAP31/CAF1 flanks) had been
  described in three papers but never deposited. SXI1 is not annotated in the
  reference assembly (W. mellicola CBS 633.66, JH668224.1); it was modelled
  with miniprot (48.7% identity to W. ichthyophaga XP_009266116.1).
- A tier-2 record built from that assembly calls all 51 BFD Wallemiales
  genomes (0 before): 33 high, 18 medium.
- SXI1, HMG and STE3 are adjacent in 33/33 genomes carrying the reference
  allele. Allele split: W. mellicola 14 vs 13; W. ichthyophaga 18 vs 4, where
  the 4 are exactly the published inverted strains.
- Update 2026-09-27: a second tier-2 putative record (W. canadensis
  EXF-10342, the other version: receptor STE3v2, no SXI1, only an HMG-box
  fragment) lets detect NAME both versions. All 51 genomes are called with a
  version: v1 33 (all high), v2 18 (all medium). Per species: W. mellicola
  14 v1 / 13 v2; W. ichthyophaga 18 / 4 (the 4 = the published inverted
  strains); W. hederae 1 / 0; W. canadensis 0 / 1. No genome carries both.
  The v2 receptor model hits the v2 genomes at 76.5-100% identity, and the
  version vote separates them widely (v2 178-234 vs v1 26-74 bits).
- Limits: no mating or meiosis has been observed in Wallemia; a near 1:1
  version split fits heterothallism but does not prove it. v1/v2 are
  placeholders; which is which mating type is unknown.
- Evidence: results/2026-09-27_puccinio_followup/, results/2026-09-27_wallemia_allele2/;
  records db/Basidiomycota/Wallemiales/671144_cbs-633-66_wallMAT_v1/ and
  1708542_exf-10342_wallMAT_v2/.

## Rhodotorula P/R and HD loci recovered genome-wide; agreement with an independent study (2026-09-27)

- Tier-2 records from the group's preprint (Liu et al., bioRxiv
  10.1101/2025.09.11.675505): 6 P/R (A1 and A2 across three clades, family
  redPR) and 5 HD (new family redHD), all on public contigs.
- 221 genomes (216 BFD Rhodotorula + 15 Sporidiobolales test genomes): HD
  0 -> 217 genomes called; P/R 215 -> 216. HD and P/R never share a contig,
  consistent with physically separate (tetrapolar) loci.
- P/R allele agrees with the study's independent strain assignments for 61 of
  62 held-out strains (66/67 including record strains); the miss is a
  reported hybrid (RIT389). A second reported hybrid is called A1+A2 with two
  HD loci.
- Limits: the comparison set is unpublished (preprint); only strains with
  public BFD assemblies were scored; the HD search costs 5.9x runtime
  (median 18 -> 106 s per genome).
- Evidence: results/2026-09-27_rhodotorula_curation/ (main checkout);
  docs/notes/2026-09-27_rhodotorula-curation.md.

## Mating pheromone receptors identified by position, from genome DNA alone (2026-09-27)

- Agaricomycete genomes carry 3-7 STE3-like receptor copies; a PF02076 tree
  does not separate mating from non-mating copies there
  (results/2026-09-27_receptor_explore/).
- Rule: a short ORF (20-130 codons, in-frame Met) ending in a strict CAAX
  motif C[VI][IV][AVMG] within 10 kb of an STE3 locus. Found from genome DNA
  with a six-frame scan; no annotation used.
- Known answers (6 genomes, 9 mating receptors, 25 other STE3 copies): strict
  CAAX flags 6/9 mating receptors and 0/25 others; chance rate 2.5% of
  random same-size windows. Every Agaricomycete mating receptor is found
  (Coprinopsis 3/3, Schizophyllum 2/2). With homology to curated pheromones:
  4/9 and 0/25, chance ~0.1%.
- Uncurated orders (52 genomes): strict CAAX flags 75 STE3 loci vs 13.6
  expected by chance (5.5x); with homology 36 vs 0.3. Polyporales 13/15 and
  Russulales 8/12 genomes flagged; flagged copies sit nearer the curated
  mating receptors in the tree (median distance 1.12-1.19 vs 2.20).
- Limits: misses Cryptococcus (pheromone 45 kb away) and Rhodotorula
  (CTIA/CTVA endings need a lineage-specific motif); fails in rusts (1/16);
  the textbook CAAX alphabet is useless (42% chance rate); ground truth is
  small (9 receptors).
- Evidence: results/2026-09-27_pheromone_positional/ (main checkout).
