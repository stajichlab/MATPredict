# MAT1-1-1 and MAT1-2-1 introns and splice classes, against the genome's own annotation (#262)

Date 2026-11-21. The census the curator asked for on 2026-10-21 (design: section 7b of `2026-10-21_gene-model-extension-design.md`), run on the Ascomycota campaign of v0.9.0. Run record and file list: `results/2026-11-21_mat_intron_census_262/RUN.txt`. Nothing in any report, call or tier changes.

## What was counted
Every MAT1-1-1 and MAT1-2-1 gene of the 19,343 reports that has an exon model (6,779 and 7,239 genes; the others, 702 and 816, are unpolished hits without a model), the introns of the model on the genome (donor and acceptor on the gene strand), and the genome's own funannotate annotation (the best-overlapping annotated mRNA: `annotated` when it covers at least half of the model, else `partial`, `novel`, or `no_annotation`). "Same intron set" means the annotated mRNA has exactly the model's introns. **This measures agreement between our model and the annotation, not which of them is right.**

## Introns per gene
| Gene | Modelled genes | 0 introns | 1 | 2 | 3 or more |
|---|---|---|---|---|---|
| MAT1-1-1 | 6,779 | 1,138 (16.8%) | 5,136 (75.8%) | 412 (6.1%) | 93 (1.4%) |
| MAT1-2-1 | 7,239 | 1,017 (14.0%) | 3,364 (46.5%) | 2,721 (37.6%) | 137 (1.9%) |

Most MAT1-1-1 genes have one intron and a third of MAT1-2-1 genes have two. The introns are short: median 52 bp (MAT1-1-1) and 53 bp (MAT1-2-1), and 90.2% (5,628 of 6,240) and 88.8% (8,201 of 9,240) are under 60 bp.

## Splice classes (per intron, on the model)
| Gene | Introns | GT-AG | GC-AG | AT-AC | other |
|---|---|---|---|---|---|
| MAT1-1-1 | 6,240 | 5,920 (94.9%) | 103 (1.7%) | 0 | 217 (3.5%) |
| MAT1-2-1 | 9,240 | 8,547 (92.5%) | 179 (1.9%) | 0 | 514 (5.6%) |

The GC-AG introns are concentrated in a few lineages: MAT1-1-1 in Pleosporales (64 of the 103; *Pyrenophora* 62 of its 76 genes) and Mycosphaerellales (23); MAT1-2-1 in Botryosphaeriales (137 of the 179). The `other` introns are the commonest donor-acceptor pairs TA-AG (36), AT-AG (19), GC-CT (18) for MAT1-1-1 and GT-TG, GT-CC, GT-AA for MAT1-2-1; most are in a few lineages (MAT1-1-1: Dothideales 40, Mycosphaerellales 31, Onygenales 28, Hypocreales 20; MAT1-2-1: Hypocreales 231, Lecanorales 97). A donor or acceptor outside the canonical set rarely occurs in nature, so these mostly mark an intron the model placed in the wrong place. Genes whose model has such an intron carry an internal stop far more often (MAT1-1-1: 32.2% against 1.9% for GT-AG-only genes).

## Agreement with the genome's annotation
| MAT1-1-1 | Genes | Annotation covers the model | Same intron set (of those annotated) | Internal stop in our model |
|---|---|---|---|---|
| no intron in our model | 1,115 | 90.6% | **6.1%** (58 of 950) | 8.6% |
| GT-AG introns only | 5,360 | 96.9% | 71.7% | 1.9% |
| has GC-AG or AT-AC | 99 | 100% | 83.5% | 2.0% |
| has a noncanonical intron | 205 | 83.1% | **4.8%** | 32.2% |

| MAT1-2-1 | Genes | Annotation covers the model | Same intron set | Internal stop |
|---|---|---|---|---|
| no intron | 979 | 88.4% | 34.8% | 11.6% |
| GT-AG only | 5,602 | 98.8% | 49.6% | 3.1% |
| has GC-AG or AT-AC | 175 | 98.8% | 3.5% | 1.7% |
| has a noncanonical intron | 483 | 87.5% | 0.3% | 6.8% |

Per intron, the exact intron is in the best-overlapping annotated mRNA: GT-AG 96.1% (MAT1-1-1) and 96.3% (MAT1-2-1); **GC-AG 82.2% of 101 (MAT1-1-1) but 5.7% of 176 (MAT1-2-1)**; noncanonical 0% of 172 and of 426. By length (MAT1-1-1 / MAT1-2-1): under 60 bp 94.0% / 93.1%, 60 to 99 bp 85.5% / 63.1%, 100 to 199 bp 73.5% / 35.6%, 200 bp or more 25.0% (20 introns) / 31.6% (19).

By lineage (`census262.by_order.tsv`, `census262.by_genus.tsv`, the n per class is in the tables):
- *Pyrenophora* MAT1-1-1: all 76 models have one intron, 62 of them GC-AG, and the annotation has exactly the same intron in 67 of 67 annotated genes. The funannotate gene set recovers the GC-AG 53 bp intron that miniprot cannot place (#262's example); the failure is miniprot's, not the annotation's.
- *Aureobasidium* MAT1-1-1: 40 of 133 introns are `other` (the TA-donor models of the issue) and the intron set matches the annotation in 46% of 139 annotated genes.
- Lineage-wide mismatches with both sides carrying an intron: *Zymoseptoria* MAT1-1-1 and Chaetothyriales MAT1-1-1 (0% same intron set), Magnaporthales and Botryosphaeriales MAT1-2-1 (0.0% and 0.2%), Glomerellales MAT1-2-1 (11.0%). Botryosphaeriales MAT1-2-1 is where the GC-AG count comes from (137 of 179), so our GC-AG models there are not in the annotation.
- Close agreement: Aspergillus (96.9% and 96.0%), Magnaporthales MAT1-1-1 (98.8%), Neurospora (98.9% and 98.7%), Parastagonospora (100%).

## Reading
- MAT1-1-1 introns are short and mostly GT-AG; the non-GT-AG share is small (5.2% of MAT1-1-1 introns, 7.5% of MAT1-2-1) but concentrated in lineages (GC-AG in Pleosporales and Botryosphaeriales, noncanonical in Dothideales, Hypocreales, Lecanorales).
- Splice class does not explain most disagreement with the annotation. A GT-AG-only model differs from the annotation in 28% (MAT1-1-1) and 50% (MAT1-2-1) of annotated genes; most disagreement is the number or position of short introns, and lineage-wide (0% in some orders, near 100% in others), which looks like a systematic difference between the two methods in those lineages. The annotation recovers 82% of MAT1-1-1 GC-AG introns, so a GC-AG intron is not by itself why a standard predictor fails; for MAT1-2-1 GC-AG the agreement is 5.7%, driven by one order.
- Intron-less models are suspect: in 94% of the MAT1-1-1 models without an intron the annotation has one, and these models carry internal stops three to five times as often as models with GT-AG introns (8.6% against 1.9%). That suggests a model that missed a short intron, or an annotation that added one that is not there; this census cannot say which.
- Models with a noncanonical intron are the least reliable (32.2% internal stops for MAT1-1-1, annotation agreement 5%), consistent with the earlier reading that they mark a misplaced intron.

## What this does not show
- Which model is right where the model and the annotation disagree: the annotation is a predictor too (not truth), and so is our exonerate or miniprot model. Truth needs RNA-seq junctions or a curated structure, as in the tools-disagree check of #349; nothing here was tested against either.
- Basidiomycota and Mucoromycotina have no MAT1-1-1 or MAT1-2-1 and were not run; the funannotate gene set is one predictor, not "gene predictors" in general (AUGUSTUS and Helixer concordance from design item 3 was not run).
- The reference distribution of design item 2 (introns of all annotated genes of the same clades, to test whether MAT genes are enriched for short or noncanonical introns) was not computed; 90% of MAT1-1-1 introns being under 60 bp is a count, not an enrichment test. The position-weight matrices for the miniprot `--spsc` stage were not built.
- Genes without an exon model (1,518) are left out of the intron counts.

Next, if wanted: an RNA-seq junction check on a sample of the disagreeing genes (the no-intron MAT1-1-1 models, and Botryosphaeriales MAT1-2-1), and the reference intron distribution.
