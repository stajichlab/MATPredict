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
