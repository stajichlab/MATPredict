# MAT intron disagreements against species-level Trinity transcripts: a pilot that mostly could not decide (#262)

Date 2026-11-21. Job 32358493 (HPCC), tables in `results/2026-11-21_mat_intron_rnaseq_262/`. Follows the intron census (#374).

## Question
Where the model and the genome's annotation disagree on introns, do RNA-seq transcripts side with the model? Two groups: MAT1-1-1 models with no intron where the annotation has one (A), and Botryosphaeriales MAT1-2-1 with a different intron set (B), plus controls where model and annotation agree (K).

## What was used
PopSeq has no RNA-seq availability for these genomes (18 RNA-Seq rows of 256,090 runs), so I used the BFD genome-guided Trinity assemblies (`Fungi_BFD_runs/rnaseq_data`), aligned with minimap2 splice mode to a window around each gene. These are species-level assemblies, so a transcript may come from a different strain. This is not the read-level junction count used for C. auris MTLA2.

## Result (287 genes)
| group | genes | Trinity file empty | fewer than 2 transcripts overlapping | with evidence |
|---|---|---|---|---|
| A, no-intron MAT1-1-1 | 212 | 167 | 202 | 10 |
| B, Botryosphaeriales MAT1-2-1 | 35 | 16 | 29 | 6 |
| K, controls | 40 | 15 | 34 | 6 |

- **The check could not decide for 265 of 287 genes.** 198 of the 287 species have a zero-byte Trinity file (no assembly), and a MAT gene is often not assembled even from a non-empty file. Of the 40 controls, 34 had too few transcripts, so no-evidence is not a sign of a bad model.
- **A, 10 genes with evidence:** in 7 (Cryomyces, Oidiodendron, Erysiphe, Opeltia, Sarocladium, Ceratocystis, Thermothelomyces) at least 2 unspliced transcripts cover the whole span and none has a junction, although the annotation has an intron. 2 are inconclusive and 1 (Oidium heveae) has a 74 bp GT-AG junction, 3 transcripts, that the model lacks. Unspliced transcripts can also come from pre-mRNA or genomic DNA, so this leans towards the model without proving it.
- **B, 6 genes with evidence:** all 6 are *Phyllosticta* (3 species, 2 genomes each, so one lineage, not six). The model's single intron has 2 transcripts in each; one of the two annotated introns has none. Two transcripts do not show the second intron is absent.
- **K, 6 genes with evidence:** 2 model_supported, 4 inconclusive on this rule. The rule is strict for agreeing genes too, so a verdict of inconclusive is weak.

## What this does and does not show
- It does not settle which model is right for the 287 genes. About 7 of 10 usable no-intron genes favour the model, which is too few to move a rule.
- For Botryosphaeriales MAT1-2-1 there is no usable evidence beyond one *Phyllosticta* lineage.
- A real answer needs reads from the sequenced strain (hisat2 on public RNA-seq, as for MTLA2), which means a per-genome ENA query. PopSeq has no such table, so the next step is either a one-off ENA exception or a PopSeq issue; that is the curator's call.
