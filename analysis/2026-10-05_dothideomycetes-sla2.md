# Dothideomycetes: SLA2 is not beside the MAT locus

Status: first test (2026-10-05), open for curator review. Prompted by the campaign dashboard: only 6% of full Dothideomycete loci carry both APN2 and SLA2.

## Question
Is the missing SLA2 flank in Dothideomycetes a detection limit (the gene is not found or modelled, or lies outside the search window) or a real difference in how the locus
is organised?

## What the v0.6.0 campaign shows
- 2,537 Dothideomycete loci (`results/2026-10-03_ascomycota_v060/loci.tsv`); all called through the generic `Ascomycota:MAT` family. Of the 2,248 full loci,
  across all 2,537 loci APN2 is present in 97% (2,468), COX13 in 89% and SLA2 in 6% (157).
- By order, loci with both flanks: Cladosporiales 63 of 91, Dothideales 34 of 215, Pleosporales 24 of 865, Mycosphaerellales 14 of 422, Botryosphaeriales 1 of 758,
  Venturiales 0 of 109.
- The database on `main` has no Dothideomycete records. Ten signed-off records (Pleosporales and Mycosphaerellales) are on the unmerged branch
  `curation-mucor-dothideo`. One of those ten mentions SLA2 and six mention APN2, so curated references would not by themselves test the question.

## Test
Genome-wide miniprot search for SLA2, APN2, COX13 and APC5 (NC1011 proteins KAI0195601.1, 603.1, 602.1, 599.1; coverage 0.5 or more, positives 0.30 or more) in 200 Dothideomycete
genomes (up to 25 per order, one called full locus each) and 60 Sordariomycete genomes as the control. For each gene the best-scoring hit is placed relative to the
called locus. Script `results/2026-10-05_dothideo_sla2_test/sla2_distance.py`; table `sla2_distance.tsv`.

## Result
| Position of the best SLA2 hit | Sordariomycetes (n=60) | Dothideomycetes (n=200) |
|---|---|---|
| inside the called locus | 80% | 14% |
| same contig, 20 to 100 kb away | 2% | 18% |
| same contig, over 100 kb away | 2% | 20% |
| another contig or no hit | 17% | 48% |
- SLA2 is found in all 200 Dothideomycete genomes, median 85% positives against the NC1011 protein (control 93%). A missed gene is not the main explanation.
- The method works where SLA2 is beside MAT: 80% inside the locus in the control.
- Within Dothideomycetes the pattern depends on the order (sampled genomes): Cladosporiales 17 of 24 inside the locus; Myriangiales 20 of 22 at 20 to 100 kb;
  Venturiales 24 of 24 on another contig; Botryosphaeriales 18 of 24 on another contig, 6 over 100 kb away; Pleosporales, Mycosphaerellales and Dothideales mostly
  on another contig or far away.
- APN2 stays with the locus (94% inside), COX13 60% (34% had no hit, a small protein that diverges), APC5 mostly absent or within 20 kb.

## Reading
Mostly real organisation, with a detection-window component. SLA2 is detached from the MAT locus in most Dothideomycete orders, by 20 to 100 kb in some
(Myriangiales) and farther or on another contig in others. The 18% at 20 to 100 kb would be found by a wider window than the 20 kb flank-carried window. Cladosporiales retain the
ancestral arrangement, so the detachment is lineage-specific, the same kind of change as in Xylariales (see `analysis/2026-10-05_xylariales-synteny.md`). Calling is not harmed: APN2
(97% of loci) serves as the flank, and the loci are still called.

## Limits
- The sample is 200 of 2,248 loci and 25 per order, so order-level shares have wide intervals (21 to 24 genomes each).
- "Another contig" includes assembly breaks and a possible wrong best hit (the control shows 17%); N50 was not checked.
- The best-scoring hit stands for the ortholog; 13 genomes have a second SLA2-like hit (the paralog family).
- Synteny is by position only; no breakpoints, and no check that the neighbouring genes are conserved around the detached SLA2.

## Next
1. Merge or review `curation-mucor-dothideo` so the campaign has Dothideomycete references (curator decision).
2. Check the genes between APN2/COX13 and the detached SLA2 in two or three orders (Myriangiales, Cladosporiales, Pleosporales) with assemblies at chromosome level.
3. Decide whether Dothideomycetes need a different flank pair in the locus model (APN2 and COX13) and whether the flank-carried window should widen for them.
