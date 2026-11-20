# The 10 long genes where the two polish tools disagree: splice-motif and model-quality check (#290, #349)

Date 2026-11-19. Item 2 of the v0.10 to-do list of #349. A recommendation per gene, not a ruling; no call, tier, label or coordinate in any report changes.

## What was checked
The 10 genes of the #290 triage (`novel_partial_proposed_actions_v090.tsv`, "curator look at the model (tools disagree)") in the v0.9.0 FungiDB run (`--preset v0.9.1`, frozen commit 7d73be5, release-71 genomes). Each carries two models in its report: the primary (`exonerate_refine` in all ten) and the `alternate_model` (`miniprot_refine`). Script `scripts/model_disagreement_check.py` (test `tests/test_model_disagreement_check.py`, 7 passed on the HPCC), job 32233314, commit 12a678d; tables in `results/2026-11-19_model_disagreement_290/` (`model_check.tsv`, one row per gene; `introns.tsv`, one row per intron of each model with its donor and acceptor; `models.faa`, both translations; `RUN.txt`). Reference lengths come from the frozen db of that run.
For each model: the introns with their donor and acceptor read on the gene's strand (`GT-AG`, `GC-AG` and `AT-AC` are canonical), the coding length, the translation (genetic code from the report: 12 for the CTG-clade yeasts), internal stops, Met start, and the length against the reference protein of the matched record. A model is CLEAN when it has no internal stop, a length divisible by 3 and only canonical introns of 30 bp or more. One clean model: that model is recommended. Both clean: both plausible. Neither: neither clean.

## Result
| Genome | Gene (tier) | Reference aa | Primary (exonerate) | Alternate (miniprot) | Recommendation |
|---|---|---|---|---|---|
| *F. odoratissimum* race 4 | MAT1-1-3 (medium) | 162 | 161 aa; 2 introns, one `GT..GG` | 182 aa; one GT-AG intron | alternate |
| *C. auris* 6684 | MTLA2 (high) | 152 | 152 aa; intron `TT..AG` | 160 aa; GT-AG | alternate |
| *C. auris* B11245 | MTLA2 (high) | 152 | 152 aa; intron `TT..AG` | 160 aa; GT-AG | alternate |
| *C. auris* B8441 | MTLA2 (high) | 152 | 152 aa; intron `TT..AG` | 160 aa; GT-AG | alternate |
| *C. pseudohaemulonii* B12108 | MTLA2 (high) | 152 | 154 aa; intron `TT..AG` | 162 aa; GT-AG | alternate |
| *H. vastatrix* race XXXIII | bE (high) | 511 | 505 aa; 3 introns, one `GT..AT` | 526 aa; 2 GT-AG introns | alternate |
| *H. vastatrix* race XXXIII | bW (high) | 494 | 316 aa, 3 internal stops; 3 introns, two noncanonical | 519 aa, 8 internal stops; one intron | neither clean |
| *P. brasiliensis* Pb18 | MAT1-1-4 (medium) | 246 | 87 aa; one GT-AG intron of 4,045 bp; no Met | 186 aa; no intron; no Met | both plausible (alternate closer in length) |
| *K. heveanensis* CBS569 | RPL39 (high) | 51 | 51 aa; 2 GT-AG introns; no Met | 58 aa; GT-AG; no Met | both plausible (primary closer) |
| *M. restricta* CBS7877 | pra1 (high) | 380 | 217 aa; GT-AG; no Met | 294 aa; GT-AG; no Met | both plausible (alternate closer; both far short of the reference) |

Totals: alternate 6, both plausible 3, neither clean 1.

## What it shows
- In all six "alternate" genes the primary model has a noncanonical intron and the miniprot model has only canonical ones. The four MTLA2 genes (three *C. auris*, *C. pseudohaemulonii*) show the same exonerate intron (`TT..AG`, 156 to 158 bp) against the same record (`498019_b11243_MTL_A`), so they look like one systematic effect of that record's alignment, not four independent errors. I did not test that.
- Where both models are clean, nothing here separates them; the length against the reference only says which is closer. The three "both plausible" genes also lack a Met start in both models, and *M. restricta* pra1 is 57% and 77% of the reference length, so both models may be truncated. *P. brasiliensis* MAT1-1-4: the primary's 4,045 bp intron is unusually long for a MAT gene (not a rule of the check), and the alternate has no intron and a lower identity (36.9% against 54.65%).
- *H. vastatrix* bW has stops and noncanonical junctions in both models and a 59 to 70% identity: it may be a disrupted copy, a poorly assembled region or a gene the reference does not fit. A curator's look at the locus, ideally with transcripts, is the only reasonable step.

## What it does not show
A canonical splice site is necessary but not sufficient: a wrong junction can sit on another GT...AG, and a real intron can be noncanonical. No RNA-seq, no protein-level comparison with a second reference and no FungiDB comparison was used (all ten are `novel`, so the annotation has no gene to compare). The check covers these ten genes only; it says nothing about how often the two tools disagree elsewhere, or which tool is right in general. Which model the report shows is unchanged.

## For v0.10, if wanted
A rule that, where the tools disagree and the primary has a noncanonical intron and the alternate has none, reports the alternate as the model would change reported gene models and coordinates in these cases. It would be a new opt-in policy (default `off`, a named preset after it passes, described in `docs/detect-options.md`), and it needs a measurement on the whole `polished_disagree` set first. Not built.
