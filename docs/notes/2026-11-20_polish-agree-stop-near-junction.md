# `polished_agree` models with a stop near a junction: a count and a direct-miniprot check (#305)

Date 2026-11-20. Measurement only; nothing in any report, call or tier changes. Recommendation at the end.

## Question
#305: the polish agreement test labels a model `polished_agree` when exonerate and miniprot agree within a tolerance on exon junctions, even when the exonerate model has an internal stop at a junction that a direct miniprot run avoids (*Wallemia ichthyophaga* STE3). How often does that happen?

## Method
Every gene of every detected locus in the three v0.9.0 definitive campaigns (`results/2026-11-09_v090_*/reports_all.tar.zst`, commit 7d73be5). Population: genes with `status: polished_agree` and `internal_stops > 0`. For each, the CDS of the stored exons was translated on the genome (genetic code from the report) and the distance from each internal stop codon to the nearest exon-exon junction measured in the CDS; "near" is a stop within 10 nt of a junction (0 = the stop codon spans it). For the near-junction genes, miniprot was run directly with the reference protein of the gene's reference record on the stored span plus or minus 5 kb (genetic code 1 for all 31), and its best overlapping model translated and checked for stops and intron motifs. Jobs 32307908 (count) and 32310115 (direct run); the two scripts were one-offs on the HPCC and are not kept. Tables and run records: `results/2026-11-20_polish_agree_stops_305/`.

## Count
| Campaign | `polished_agree` genes | with an internal stop | stop within 10 nt of a junction |
|---|---|---|---|
| Ascomycota | 29,169 | 128 | 10 (plus 1 not checked: MTLA1 of GCA_003706265.2, genetic code 26, which the check does not support) |
| Basidiomycota | 3,399 | 23 | 19 |
| Mucoromycotina | 1,806 | 6 | 2 |
| All | 34,374 | 157 (0.46%) | 31 (0.09%) |

Of the 31, the stop spans a junction (distance 0) in 24, and lies 1 to 7 nt from one in 7. By gene: *Wallemia* STE3 17 (all 17 genomes of that species, 18 of the 19 Basidiomycota rows are at distance 0), SLA2 5 (the MATtub family), MAT1-2-1 4, and one each of MAT1-2-10, bW, beta_fg, rnhA and algA. For comparison, `polished_disagree` has 3,482 genes with an internal stop and `polished_single` 455 (`*.status.tsv`).

## Direct miniprot on the 31
| Result of the direct miniprot model | Genes |
|---|---|
| stop-free, canonical introns | 24: *Wallemia* STE3 17 (360 aa, 3 GT-AG introns, against the stored 357 aa with one stop), bW, beta_fg, rnhA, algA and 3 MAT1-2-1 |
| one internal stop left | 2 (MAT1-2-10 and one MAT1-2-1) |
| miniprot reports a frameshift | 5 SLA2 (MATtub genomes) |

The five SLA2 models translate with 66 to 70 stops because miniprot reports `Frameshift=1` for the alignment (checked on GCA_024713625.1: the 5' exon aligns at 85% identity with a one-base gap, and the other four have the same stop signature); my check concatenates the exons without applying the frameshift, so these are not stop-free models and not comparable. I did not look further into them. SLA2 is a flank gene, not a MAT gene, so no call depends on it.

## What this says
- The #305 failure exists, is rare (0.09% of `polished_agree` genes), and in these runs is mostly one lineage: *Wallemia* STE3 is 17 of 31, one species and one allele, so these are not 17 independent observations. In 24 of 31 a direct miniprot run gives a stop-free model with canonical introns.
- Nothing here shows that the tolerance should shrink: I did not measure how many correct `polished_agree` models a smaller tolerance would flip, and the rare stops are found by the stop count anyway.
- The 157 genes are already marked `defective` in the reports (internal stop, issue #151), so a reader is told the model has a stop. Whether the stop is real (a pseudogene) or an artefact of the model was not tested; no RNA-seq or genome-annotation check was made.
- This census does not cover genes without a stop where the models agree but are both wrong, or `polished_single` genes.

## Recommendation (not a ruling)
1. Do not shrink the agreement tolerance. It is untested, and the targeted case is 0.09% of the label.
2. If the 24 genes are worth fixing, add an opt-in re-test: when a `polished_agree` model has an internal stop within 10 nt of a junction, run miniprot directly and take its model only if it is stop-free with canonical introns, and label the gene `junction_artifact`. Default off, in no preset, described in `docs/detect-options.md`, with the usual before/after on the 510 panel genomes and the 313 FungiDB genomes (CLAUDE.md). I would not build it for 31 genes in 34,374 unless the *Wallemia* STE3 case matters to a call or a table.
3. For *Wallemia* STE3 itself the cheaper route may be the reference record: the stored model matches the reference at about 65% identity and the direct model is 3 aa longer than the stored one (360 against 357) in every genome. A curator can look at `671144_cbs-633-66_wallMAT_v1` (STE3 366 aa in the record).
