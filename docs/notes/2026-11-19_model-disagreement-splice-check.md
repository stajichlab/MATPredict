# Genes where the two polish tools disagree: splice-motif and model-quality check, and a census of all of them (#290, #349)

Date 2026-11-19. Item 2 of the v0.10 to-do list of #349. Recommendations only; no call, tier, label or coordinate in any report changes.

**Correction.** The first version of this note (the same day) recommended the miniprot model for six of the ten genes on splice motifs alone. The census below shows that the motifs and the reference disagree for four of them, so those are now reported as a conflict. That earlier recommendation is retracted.

## 1. The 10 genes of the #290 triage
The ten genes of `novel_partial_proposed_actions_v090.tsv` marked "tools disagree" in the v0.9.0 FungiDB run (`--preset v0.9.1`, frozen commit 7d73be5). Each report carries the primary model (`exonerate_refine` in all ten) and the `alternate_model` (`miniprot_refine`). `scripts/model_disagreement_check.py` (test `tests/test_model_disagreement_check.py`) reads the introns of both models on the genome with the donor and acceptor dinucleotides on the gene's strand (`GT-AG`, `GC-AG`, `AT-AC` canonical), translates them (genetic code from the report), and counts internal stops. A model is CLEAN with no internal stop, a length divisible by 3 and only canonical introns of 30 bp or more. It then recommends: one clean model, that model; both clean, both plausible; neither, neither clean; and `conflict` when the model the motifs favour has an identity more than 5 points lower than the other AND is farther from the reference protein's length. Job 32236188, commit 7130123; tables `results/2026-11-19_model_disagreement_290/{model_check,introns}.tsv`, `models.faa`.

| Genome | Gene (tier) | Ref aa | Primary (exonerate): aa, identity, introns | Alternate (miniprot): aa, identity, introns | Recommendation |
|---|---|---|---|---|---|
| *F. odoratissimum* race 4 | MAT1-1-3 (medium) | 162 | 161, 88.75, one `GT..GG` + GT-AG | 182, 79.12, GT-AG | conflict |
| *C. auris* 6684 | MTLA2 (high) | 152 | 152, 98.68, `TT..AG` | 160, 92.5, GT-AG | conflict |
| *C. auris* B11245 | MTLA2 (high) | 152 | 152, 100.0, `TT..AG` | 160, 93.75, GT-AG | conflict |
| *C. auris* B8441 | MTLA2 (high) | 152 | 152, 98.68, `TT..AG` | 160, 92.5, GT-AG | conflict |
| *C. pseudohaemulonii* B12108 | MTLA2 (high) | 152 | 154, 52.0, `TT..AG` | 162, 47.24, GT-AG | alternate |
| *H. vastatrix* race XXXIII | bE (high) | 511 | 505, 63.17, one `GT..AT` + 2 GT-AG | 526, 60.08, 2 GT-AG | alternate |
| *H. vastatrix* race XXXIII | bW (high) | 494 | 316, 69.87, 3 stops, 2 noncanonical | 519, 59.36, 8 stops | neither clean |
| *P. brasiliensis* Pb18 | MAT1-1-4 (medium) | 246 | 87, 54.65, GT-AG of 4,045 bp; no Met | 186, 36.9, no intron; no Met | both plausible |
| *K. heveanensis* CBS569 | RPL39 (high) | 51 | 51, 90.0; no Met | 58, 74.14; no Met | both plausible |
| *M. restricta* CBS7877 | pra1 (high) | 380 | 217, 24.07; no Met | 294, 23.76; no Met | both plausible |

Totals: conflict 4, alternate 2, both plausible 3, neither clean 1. In the four conflicts the exonerate model reproduces the reference at 98.7 to 100% identity and (for the *C. auris* MTLA2 genes) at exactly the reference length, with an unusual `TT..AG` intron of 158 bp. That is either a real noncanonical intron or a reference-derived model; these checks cannot tell. The three "both plausible" genes have no Met start in either model and *M. restricta* pra1 is far shorter than its reference in both. *H. vastatrix* bW has stops and noncanonical junctions in both and needs a curator's look at the locus.

## 2. Census of every `polished_disagree` gene (v0.9.0 definitive run)
`scripts/model_disagreement_census.py` (test `tests/test_model_disagreement_census.py`) applies the same CLEAN test to both models of every `polished_disagree` gene of every detected locus in the three campaigns (`reports_all.tar.zst`, genomes from `input_clean_genomes/clean`, or the campaign `inputs.tsv` for the Mucoromycotina samples named `<source>__<id>`). Job array 32233690 (Ascomycota, Basidiomycota) and 32235022 (Mucoromycotina, rerun with the genome-path fix after a first run could not find the prefixed samples); tables `census/<run>.summary.tsv` (the per-gene `census.tsv` files, 48,000 rows, stay on the HPCC in `results/2026-11-19_model_disagreement_290/census/`).

| Run | Genes checked | Both clean | Primary issue, alternate clean (the rule's switches) | Primary clean, alternate issue | Both issue |
|---|---|---|---|---|---|
| Ascomycota | 37,771 (18 not checked) | 21,302 | 9,727 (26%) | 2,978 | 3,764 |
| Basidiomycota | 8,748 | 4,756 | 1,832 (21%) | 1,368 | 792 |
| Mucoromycotina | 1,539 | 787 | 422 (27%) | 30 | 300 |
| **All** | **48,058** | **26,845** | **11,981 (25%)** | **4,376** | **4,856** |

A rule "prefer the alternate when the primary has a noncanonical intron and the alternate has none" would act on 10,534 + 2,245 + 457 = 13,236 genes (the trigger counts); the switches above are the subset where the alternate is also clean. For the 11,981 switches, the gene span would change in 9,201 + 1,745 + 353 = 11,299 and the alternate lies partly outside the reported locus span in 2,742 + 516 + 50 = 3,308.

**Yardsticks for the switches** (`scripts/model_disagreement_summary.py`, `yardsticks.tsv`; identity and length come from the tools that made the models, so they are not independent of them, but they are independent of the splice-motif test; exonerate aligns against the reference and favours its own model, so the length column is the less biased one):

| Run | Switch genes | Alternate has higher identity | Median identity difference (alt minus primary) | Alternate closer to the reference length |
|---|---|---|---|---|
| Ascomycota | 9,727 | 2,126 (22%) | -2.7 points | 3,524 of 9,720 (36%) |
| Ascomycota, core MAT genes | 2,477 | 262 (11%) | -6.2 | 582 of 2,470 (24%) |
| Basidiomycota | 1,832 | 240 (13%) | -3.0 | 959 of 1,794 (53%) |
| Basidiomycota, core MAT genes | 891 | 127 (14%) | -5.3 | 594 of 853 (70%) |
| Mucoromycotina | 422 | 111 (26%) | -0.9 | 202 of 422 (48%) |

## What this shows
- The tools disagree on a quarter of a million... no: on 48,058 genes, and in one quarter of those only the primary has a splice-motif issue. The reverse (primary clean, alternate with an issue) is common too (4,376).
- A noncanonical intron in the primary is not a reliable sign that the primary is the worse model. In the largest group, Ascomycota, the alternate is lower in identity in 78% of the switches and farther from the reference length in 59%; in Basidiomycota core genes the alternate is closer in length in 70%, so the picture differs by clade. A single rule over all clades would trade one set of errors for another.
- The FungiDB annotation cannot settle it: of the 87 `polished_disagree` genes in the 313 FungiDB genomes (`annotation/fungidb_v090.summary.tsv`, `scripts/model_vs_annotation.py`), 66 are in genomes with no GFF, 20 have no overlapping annotated gene, and 1 overlaps and matches neither model.

## What it does not show
Which model is right in any single gene. A canonical splice site is necessary, not sufficient. No RNA-seq, no second reference and no third aligner were used. The census says nothing about genomes outside these campaigns.

## Recommendation for v0.10 (not a ruling)
Do not build "prefer the alternate when the primary has a noncanonical intron" as stated: it would change the gene span of about 11,300 genes in a v0.9-sized run, and by the yardsticks above it would more often move to the model that is farther from the reference. What would justify any such rule is independent evidence on a sample, for example splice junctions from RNA-seq on genomes that have public RNA-seq, or a third aligner, scored over a few hundred disagreeing genes stratified by clade. If a narrower rule is still wanted, it would need to require that the alternate also has equal or higher identity and is not farther from the reference length, and it must be an opt-in policy (default `off`, described in `docs/detect-options.md`, with the acceptance test of `CLAUDE.md`). The ten genes of section 1 do not need a rule: the four conflicts and *H. vastatrix* bW are for a curator's look.
