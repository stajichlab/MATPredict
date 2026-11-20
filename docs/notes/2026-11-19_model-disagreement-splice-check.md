# Genes where the two polish tools disagree on the model: ten genes for curation, and why there is no rule (#290, #349)

Date 2026-11-19. Item 2 of the v0.10 list of #349, closed as curation only (curator, 2026-11-19). Nothing in any report, call or tier changes.

## What "disagree" means here
`polished_disagree` means both `exonerate --refine` and `miniprot` modelled the gene and the exon count or a boundary differs by more than 10 bp. Exonerate's model is canonical by design (`polish.classify`); miniprot's is kept as `alternate_model`. The status does not affect tiering (`pipeline.py`: "both are polished as far as tiering is concerned"), a call or an idiomorph; it affects only the exon/intron structure and protein reported for that gene. The one exception is the opt-in weak-fragment cap, which ranks "agree" above "disagree" (in no preset).

## The ten genes of the #290 triage
Genes marked "tools disagree" in `novel_partial_proposed_actions_v090.tsv`, v0.9.0 FungiDB run (`--preset v0.9.1`, commit 7d73be5). For each model: introns with donor and acceptor on the gene's strand (canonical `GT-AG`, `GC-AG`, `AT-AC`), internal stops, identity to the reference and length against the reference protein. Tables: `results/2026-11-19_model_disagreement_290/{model_check,introns}.tsv`, `models.faa`, `RUN.txt`. The code that made them was a one-off, run on the HPCC at commit 7130123 of the PR #363 branch history, and is not kept in the repository.

| Genome | Gene (tier) | Ref aa | Primary (exonerate): aa, identity, introns | Alternate (miniprot): aa, identity, introns | Reading |
|---|---|---|---|---|---|
| *F. odoratissimum* race 4 | MAT1-1-3 (medium) | 162 | 161, 88.75, one `GT..GG` + GT-AG | 182, 79.12, GT-AG | conflict; RNA-seq uninformative |
| *C. auris* 6684 | MTLA2 (high) | 152 | 152, 98.68, `TT..AG` | 160, 92.5, GT-AG | conflict; RNA-seq supports miniprot |
| *C. auris* B11245 | MTLA2 (high) | 152 | 152, 100.0, `TT..AG` | 160, 93.75, GT-AG | conflict; RNA-seq (other strains) supports miniprot |
| *C. auris* B8441 | MTLA2 (high) | 152 | 152, 98.68, `TT..AG` | 160, 92.5, GT-AG | conflict; RNA-seq supports miniprot |
| *C. pseudohaemulonii* B12108 | MTLA2 (high) | 152 | 154, 52.0, `TT..AG` | 162, 47.24, GT-AG | motifs favour miniprot; identities close |
| *H. vastatrix* race XXXIII | bE (high) | 511 | 505, 63.17, one `GT..AT` + 2 GT-AG | 526, 60.08, 2 GT-AG | motifs favour miniprot; identities close |
| *H. vastatrix* race XXXIII | bW (high) | 494 | 316, 69.87, 3 stops, 2 noncanonical | 519, 59.36, 8 stops | neither clean |
| *P. brasiliensis* Pb18 | MAT1-1-4 (medium) | 246 | 87, 54.65, GT-AG of 4,045 bp; no Met | 186, 36.9, no intron; no Met | both plausible |
| *K. heveanensis* CBS569 | RPL39 (high) | 51 | 51, 90.0; no Met | 58, 74.14; no Met | both plausible |
| *M. restricta* CBS7877 | pra1 (high) | 380 | 217, 24.07; no Met | 294, 23.76; no Met | both plausible |

"Conflict" = the splice motifs favour miniprot (the primary has a noncanonical intron) but the primary has more than 5 points higher identity and is closer to the reference length. In the four conflicts the exonerate model reproduces the reference at 98.7 to 100% identity (for the *C. auris* MTLA2 genes at exactly the reference length) with an unusual `TT..AG` intron of 158 bp. Splice motifs, identity and length are not independent evidence of which model is right; the next section is.

## RNA-seq junctions for the four conflicts
**Where public RNA-seq exists** (one-off read-only ENA lookup, 2026-11-19, `library_source="TRANSCRIPTOMIC"`; PopSeq holds only WGS/GENOMIC runs, so this is an exception to the PopSeq rule, approved by the curator for this check; `rnaseq_oneoff/selected_studies_runs.tsv`): *C. auris* 1,143 runs, with the exact strains B8441 (PRJNA445471, PRJNA682185, PRJNA1330823; 63 runs) and 6684 (strain Ci6684, PRJNA510662; 32 runs) and none for B11245; *F. odoratissimum* 26 runs, among them the TR4 isolates B2 (PRJNA449165, PRJNA478712) and Foc028 (PRJNA1273293). None is a sexual-stage or mating-type study.

**Method** (one-off job array 32237138, rerun of the *F. odoratissimum* task as 32239288 after a single-end file name was not expected; the script is `rnaseq_junctions_oneoff/rnajunc_oneoff.sh` on the HPCC, not in the repository): three untreated or control runs per genome (two for B11245, which has no RNA-seq of its own, so B8441 and 6684 reads were mapped to it: other-strain evidence), read counts checked against ENA, `hisat2 --dta --max-intronlen 5000`, then for each model's introns the number of alignments whose CIGAR skips exactly that intron, and every other junction with at least 5 reads in the region. Alignments are counted, not read pairs. Tables: `rnaseq_oneoff/junctions_<sample>.tsv`, `RUN_<sample>.txt`.

| Genome (runs) | Alternate (miniprot) intron | Reads across it, per run | Primary (exonerate) intron | Reads across it, per run |
|---|---|---|---|---|
| *C. auris* B8441 (3 controls, 55 to 60% of reads align) | 852688-852821, 134 bp, GT-AG | 10, 8, 4 | 852693-852850, 158 bp, `TT..AG` | 0, 0, 0 |
| *C. auris* 6684 (3 wild-type, 98% align) | 179576-179709, 134 bp, GT-AG | 3, 11, 6 | 179581-179738, 158 bp, `TT..AG` | 0, 0, 0 |
| *C. auris* B11245 (B8441 and 6684 reads, other strains) | 1585457-1585590, 134 bp, GT-AG | 10, 1 | 1585462-1585619, 158 bp, `TT..AG` | 0, 0 |
| *F. odoratissimum* race 4 (B2, B2, Foc028; about 169M reads) | shared with the primary: 339814-339860, 47 bp, GT-AG | 0, 1, 0 | primary only: 339561-339620, 60 bp, `GT..GG` | 0, 0, 0 |

**Reading.** For the three *C. auris* MTLA2 genes the RNA-seq supports the miniprot model: every run has reads across its 134 bp GT-AG intron and none across the exonerate 158 bp `TT..AG` intron, and no other junction reaches 5 reads. Counts are low (MTLA2 is lowly expressed), the B8441 libraries map poorly (55 to 60%), and B11245 rests on other-strain reads, so this is strong for the direction and weak in depth. For *F. odoratissimum* MAT1-1-3 the RNA-seq is uninformative: the gene is essentially unexpressed in these conditions (one read on the shared intron in about 169M reads), so it neither supports nor rules out the noncanonical primary intron. The other six genes were not checked for RNA-seq (the taxon identifiers I first tried for *C. pseudohaemulonii* and *H. vastatrix* were guesses and are not used).

**Why the identity and length yardsticks failed here.** The reference record of these genes, `db/Ascomycota/Serinales/498019_b11243_MTL_A`, models MTLA2 with exons 14845-15087 and 15246-15461: the depositors' annotation of B11243 (PSK75931.1, 152 aa, tier 2), whose intron is 158 bp, the intron the exonerate model has. Exonerate aligns to that record, so its model reproduces the record at 100% identity and 152 aa by construction; the RNA-seq says the record's intron is probably wrong. A reference-based yardstick is only as good as the reference's gene model, and it is circular where the reference carries the same error as the primary. A correction to the record is a curation proposal (a curator accepts `db/` changes) and is listed below, not made.

## Census of all `polished_disagree` genes (v0.9.0 definitive run)
The same test applied to both models of every `polished_disagree` gene of the three campaigns (job arrays 32233690 and 32235022, HPCC; per-gene tables not kept). CLEAN = no internal stop, length divisible by 3, canonical introns of 30 bp or more.

| Run | Genes | Both clean | Primary issue, alternate clean | Primary clean, alternate issue | Both issue |
|---|---|---|---|---|---|
| Ascomycota | 37,771 (18 not checked) | 21,302 | 9,727 (26%) | 2,978 | 3,764 |
| Basidiomycota | 8,748 | 4,756 | 1,832 (21%) | 1,368 | 792 |
| Mucoromycotina | 1,539 | 787 | 422 (27%) | 30 | 300 |
| All | 48,058 | 26,845 | 11,981 (25%) | 4,376 | 4,856 |

For the 11,981 genes where only the primary has an issue ("switches"), the alternate has the higher identity to the reference in 2,126 of 9,727 (22%) in Ascomycota, 240 of 1,832 (13%) in Basidiomycota and 111 of 422 (26%) in Mucoromycotina, and is closer to the reference length in 36%, 53% and 48% (Ascomycota core MAT genes 24%, Basidiomycota core MAT genes 70%). Switching would change the gene span of 11,299 genes and put the alternate partly outside the reported locus in 3,308. Identity favours exonerate (it aligns against the reference), and, as the *C. auris* MTLA2 case shows, a reference that carries the same error as the primary makes both identity and length circular, so these two yardsticks are weak evidence in either direction. The FungiDB annotation cannot decide it: of 87 such genes in the 313 FungiDB genomes, 66 are in genomes with no GFF, 20 have no overlapping annotated gene and 1 overlaps and matches neither model.

## Decision (revised after the RNA-seq)
Item 2 stays closed as curation only (curator, 2026-11-19): no pipeline rule. What changes is the reason. The census shows the tools disagree on tens of thousands of genes and that a motif-based switch would change about 11,300 gene spans, but its yardsticks (identity, length) are circular where the reference carries the primary's error, so they do not show that such a rule picks the worse model. The only independent evidence so far (RNA-seq, three genomes of one gene) favours the motif-clean miniprot model. A policy would need a stratified sample of disagreeing genes scored against RNA-seq junctions, which is a larger task not started here.

**For a curator**
- *C. auris* MTLA2: propose a correction of the exon structure of `498019_b11243_MTL_A` (MTLA2: the 134 bp GT-AG intron, protein about 160 aa, in place of the 158 bp `TT..AG` intron, 152 aa). Whether to accept it, and whether the other *Candidozyma*/Serinales MTLA2 records have the same intron, are the curator's calls. It would probably also remove the `polished_disagree` status of these three genes in a later run; that was not tested.
- *F. odoratissimum* MAT1-1-3: needs RNA-seq from a condition where the gene is expressed, or a sequence-based argument; nothing decided.
- *H. vastatrix* bW (neither model clean), *P. brasiliensis* MAT1-1-4, *K. heveanensis* RPL39, *M. restricta* pra1 (both plausible, no Met start): the check gives no answer; RNA-seq was not looked up for them.
- *C. pseudohaemulonii* MTLA2 and *H. vastatrix* bE (motifs favour miniprot; identities close): no RNA-seq checked.
