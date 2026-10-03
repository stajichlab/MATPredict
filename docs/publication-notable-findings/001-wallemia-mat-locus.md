# 001. A putative Wallemia MAT locus detected in all 51 genomes, with two named alleles

- **Category:** curation-first; biology
- **Status:** candidate (the locus is putative: no mating has been observed in *Wallemia*)
- **Lineage:** Basidiomycota, Wallemiomycotina, Wallemiales

## Summary
A putative MAT locus (SXI1, HMG, STE3, with BAP31 and CAF1 flanks) was
described in three papers but never deposited, and SXI1 is not annotated in the
reference assembly. Two tier-2 records, one per allele, call all 51 BFD
Wallemiales genomes (0 before) and name the allele in each. The allele split
reproduces the published inverted strains exactly.

## Evidence
- Record `671144_cbs-633-66_wallMAT_v1` (allele v1), *W. mellicola* CBS 633.66,
  `JH668224.1:17,113-42,263`: BAP31 EIM23475.1, STE3 EIM23479.1, CAF1
  EIM23480.1, HMG EIM23483.1 (all 100% to their accessions); SXI1 unannotated,
  modelled at 41,384-41,863 + 41,907-42,263 (miniprot from *W. ichthyophaga*
  XP_009266116.1, 48.7% over 263 aa; 278 aa). Source:
  `results/2026-09-27_puccinio_followup/NOTE.md` (main checkout).
- Record `1708542_exf-10342_wallMAT_v2` (allele v2), *W. canadensis*
  EXF-10342, GCA_056320075.1, `JBHFOR010000004.1:328,288-351,767`: STE3v2
  KAN3001147.1, CAF1 KAN3001148.1, BAP31 KAN3001158.1 (100%); HMG and SXI1
  absent (43-aa HMG-box fragment only). Source:
  `results/2026-09-27_wallemia_allele2/NOTE.md` (main checkout).
- Detection: 51/51 called; 33 v1 (high), 18 v2 (medium before the
  2026-09-27 tier change); no genome gets both alleles. v2 best-bitscore
  178-234 vs v1 26-74 in v2 genomes.
- Synteny (independent of detection): SXI1, HMG and STE3 on one contig within
  50 kb in 33/33 v1 genomes; all five genes within 60 kb in 32/33.
- Allele split: *W. mellicola* 14 v1 / 13 v2; *W. ichthyophaga* 18 / 4, the
  4 = EXF-759, EXF-3555, EXF-8622, EXF-8623 (the inverted strains of
  Gostincar et al. 2019); *W. hederae* 1 / 0; *W. canadensis* 0 / 1.
- Within-species STE3 identity: v1 vs v2 about 29-33%, v1 vs v1 65-100%.

## Method that found it
Tier-2 records built from genome annotation plus miniprot (curation-puccinio
commits `6590929`, `3d8a755`); detection with a per-allele receptor split in
family `wallMAT` (enum v1/v2); synteny and allele checks by tblastn and
miniprot (`checks.py`, `checks_miniprot.py`, `versions.tsv` in the follow-up
folder). Curator sign-off 2026-09-27 (`4d454f3`, `ac33880`).

## Verification done / still open
- Done: synteny in 33/33; allele split matches Gostincar et al. 2019 and the
  about-half split of Sun et al. 2019 (per the follow-up note).
- Open: which allele is which mating type; any evidence of mating or meiosis;
  SXI1 HD1/HD2 class (deferred until a tree).

## Limits
- The locus is putative in the literature. A 1:1 split fits heterothallism
  but does not prove it.
- SXI1 is this project's model, not a deposited gene. STE3v2 (192 aa) is
  probably 3'-truncated (annotation report B5.5).

## Related
Annotation report B5.4, B5.5; `docs/notes/2026-09-26_pucciniomycotina-curation.md` (curation-puccinio).
