# Genes where the two polish tools disagree on the model: ten genes for curation, and why there is no rule (#290, #349)

Date 2026-11-19. Item 2 of the v0.10 list of #349, closed as curation only (curator, 2026-11-19). Nothing in any report, call or tier changes.

## What "disagree" means here
`polished_disagree` means both `exonerate --refine` and `miniprot` modelled the gene and the exon count or a boundary differs by more than 10 bp. Exonerate's model is canonical by design (`polish.classify`); miniprot's is kept as `alternate_model`. The status does not affect tiering (`pipeline.py`: "both are polished as far as tiering is concerned"), a call or an idiomorph; it affects only the exon/intron structure and protein reported for that gene. The one exception is the opt-in weak-fragment cap, which ranks "agree" above "disagree" (in no preset).

## The ten genes of the #290 triage
Genes marked "tools disagree" in `novel_partial_proposed_actions_v090.tsv`, v0.9.0 FungiDB run (`--preset v0.9.1`, commit 7d73be5). For each model: introns with donor and acceptor on the gene's strand (canonical `GT-AG`, `GC-AG`, `AT-AC`), internal stops, identity to the reference and length against the reference protein. Tables: `results/2026-11-19_model_disagreement_290/{model_check,introns}.tsv`, `models.faa`, `RUN.txt`. The code that made them was a one-off, run on the HPCC at commit 7130123 of the PR #363 branch history, and is not kept in the repository.

| Genome | Gene (tier) | Ref aa | Primary (exonerate): aa, identity, introns | Alternate (miniprot): aa, identity, introns | Reading |
|---|---|---|---|---|---|
| *F. odoratissimum* race 4 | MAT1-1-3 (medium) | 162 | 161, 88.75, one `GT..GG` + GT-AG | 182, 79.12, GT-AG | conflict |
| *C. auris* 6684 | MTLA2 (high) | 152 | 152, 98.68, `TT..AG` | 160, 92.5, GT-AG | conflict |
| *C. auris* B11245 | MTLA2 (high) | 152 | 152, 100.0, `TT..AG` | 160, 93.75, GT-AG | conflict |
| *C. auris* B8441 | MTLA2 (high) | 152 | 152, 98.68, `TT..AG` | 160, 92.5, GT-AG | conflict |
| *C. pseudohaemulonii* B12108 | MTLA2 (high) | 152 | 154, 52.0, `TT..AG` | 162, 47.24, GT-AG | motifs favour miniprot; identities close |
| *H. vastatrix* race XXXIII | bE (high) | 511 | 505, 63.17, one `GT..AT` + 2 GT-AG | 526, 60.08, 2 GT-AG | motifs favour miniprot; identities close |
| *H. vastatrix* race XXXIII | bW (high) | 494 | 316, 69.87, 3 stops, 2 noncanonical | 519, 59.36, 8 stops | neither clean |
| *P. brasiliensis* Pb18 | MAT1-1-4 (medium) | 246 | 87, 54.65, GT-AG of 4,045 bp; no Met | 186, 36.9, no intron; no Met | both plausible |
| *K. heveanensis* CBS569 | RPL39 (high) | 51 | 51, 90.0; no Met | 58, 74.14; no Met | both plausible |
| *M. restricta* CBS7877 | pra1 (high) | 380 | 217, 24.07; no Met | 294, 23.76; no Met | both plausible |

"Conflict" = the splice motifs favour miniprot (the primary has a noncanonical intron) but the primary has more than 5 points higher identity and is closer to the reference length. In the four conflicts the exonerate model reproduces the reference at 98.7 to 100% identity (for the *C. auris* MTLA2 genes at exactly the reference length) with an unusual `TT..AG` intron of 158 bp: a real noncanonical intron or a model that follows the reference; these checks cannot tell, and an earlier recommendation of the miniprot model for them, made on motifs alone, is retracted. Splice motifs, identity and length are not independent evidence of which model is right. Real evidence would be splice junctions from RNA-seq (see the RNA-seq section, if added below).

## Census of all `polished_disagree` genes (v0.9.0 definitive run)
The same test applied to both models of every `polished_disagree` gene of the three campaigns (job arrays 32233690 and 32235022, HPCC; per-gene tables not kept). CLEAN = no internal stop, length divisible by 3, canonical introns of 30 bp or more.

| Run | Genes | Both clean | Primary issue, alternate clean | Primary clean, alternate issue | Both issue |
|---|---|---|---|---|---|
| Ascomycota | 37,771 (18 not checked) | 21,302 | 9,727 (26%) | 2,978 | 3,764 |
| Basidiomycota | 8,748 | 4,756 | 1,832 (21%) | 1,368 | 792 |
| Mucoromycotina | 1,539 | 787 | 422 (27%) | 30 | 300 |
| All | 48,058 | 26,845 | 11,981 (25%) | 4,376 | 4,856 |

For the 11,981 genes where only the primary has an issue ("switches"), the alternate has the higher identity to the reference in 2,126 of 9,727 (22%) in Ascomycota, 240 of 1,832 (13%) in Basidiomycota and 111 of 422 (26%) in Mucoromycotina, and is closer to the reference length in 36%, 53% and 48% (Ascomycota core MAT genes 24%, Basidiomycota core MAT genes 70%). Switching would change the gene span of 11,299 genes and put the alternate partly outside the reported locus in 3,308. Identity favours exonerate (it aligns against the reference), so length is the less biased of the two. The FungiDB annotation cannot decide it: of 87 such genes in the 313 FungiDB genomes, 66 are in genomes with no GFF, 20 have no overlapping annotated gene and 1 overlaps and matches neither model.

## Decision
No rule. A "prefer the alternate when the primary has a noncanonical intron" policy would change about 11,300 gene spans and, in the largest clade, more often move to the model that is lower in identity and farther from the reference. The choice affects only the reported gene structure, not calls or tiers. The ten genes are for a curator's look, with RNA-seq where it exists.
