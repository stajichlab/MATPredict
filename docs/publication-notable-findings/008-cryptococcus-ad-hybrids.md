# 008. Cryptococcus AD hybrids called a+alpha

- **Category:** hybrid
- **Status:** verified for CBS 132 (known AD hybrid, NCBI organism name); curator identification for PH4021C
- **Lineage:** Basidiomycota, Tremellomycetes, Tremellales

## Summary
With the SXI1/SXI2 bonus slot, two C. neoformans genomes change from alpha to
a+alpha. Both are AD hybrids, so two idiomorphs are expected and are not
evidence of homothallism.

## Evidence
- CBS 132, GCA_006992865.1 (NCBI organism "Cryptococcus neoformans x
  Cryptococcus deneoformans"; 27.2 Mb, 559 contigs, Canu, PRJNA477598).
- PH4021C, GCA_025531435.1 (28.0 Mb, 530 contigs, SPAdes, PRJNA833221; a
  haploid C. neoformans is about 19 Mb). Calls: alpha high `mat_locus` on
  JAMYHK010000168.1 (FAO1, STE3, SXI1); a medium `partial_locus` on
  JAMYHK010000139.1 (FAO1, SXI2, RPL39).
- Genotype changes (alpha -> a+alpha) for exactly these two:
  `results/2026-09-26_sxi_slot_gateA/compare_9655f42.txt` (committed).
- Assembly sizes: NCBI Datasets dataset_report, queried 2026-09-27.

## Method that found it
SXI1/SXI2 bonus slot (commit `9655f42`), 243-genome Cryptococcus re-run.

## Verification done / still open
- CBS 132: hybrid per NCBI taxonomy. PH4021C: identified as an AD hybrid by
  the curator (2026-09-26); the assembly size fits but no subgenome test was run.
- Open: a subgenome check for PH4021C (A vs D identity of the a-locus genes).

## Limits
PH4021C's hybrid status rests on curator knowledge and assembly size.

## Related
Entry 010 (Rhodotorula hybrids).
