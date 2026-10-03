# Study index

Status: decided = ruling made; open = decision pending; running = results not
in yet; superseded = replaced by a later study. `results/` = main checkout.

## Mucoromycota
| Date | Report | Status | Key numbers |
|---|---|---|---|
| 09-26 | [Early-diverging scan and synteny](2026-09-26_early-diverging-scan-and-synteny.md) | decided | Mucoromycota 227/293; Mortierello/Kickxello 0 synteny support |
| 09-26 | [sexM/sexP trees](2026-09-26_sexMP-trees.md) | decided | label/clade 204/206 (FastTree); final tree sexP UFBoot 62 |
| 09-26 | [Idiomorph labelling and classifier](2026-09-26_idiomorph-labelling-and-classifier.md) | decided | btbA artefact; LOO 85/85; Zygo 23/23 |
| 09-26 | [Relaxed pass / Phycomyces](2026-09-26_relaxed-pass-phycomyces.md) | decided | defect since 6bc985d fixed |
| 09-27 | [Umbelopsis, guard, Circinella](2026-09-27_umbelopsis.md) | decided; merged a2fe1b4 | 13/14 called; records signed off; guard dropped; Circinella lost to classifier rebuild |
| 09-27 | [R. arrhizus split locus](2026-09-27_rarrhizus-split-locus.md) | decided | 11 of 12 GL genomes recovered (corrected from 12) |
| 09-29 | [Classifier builds: aligner, determinism, gate, P1](2026-09-29_classifier-builds.md) | decided; shipped 9c39e2c | no aligner more accurate; mafft --auto = L-INS-i; ClipKIT 16 wrong; builds byte-identical; P1 withholds 26, reveals 5 |
| 09-29 | [sexM-like paralogs, strain labels, Absidia](2026-09-29_paralogs-and-strain-labels.md) | decided | labels 10/7/4 (n=21); P1 in both mating types; Absidia flanks off-scaffold |
| 09-29 | [Two idiomorphs and homothallism](2026-09-29_two-idiomorphs-and-homothallism.md) | decided (report-only); causes open | 36 unlinked; literature signal 13/16 vs 5/9; Syzygites both idiomorphs |

## Basidiomycota
| Date | Report | Status | Key numbers |
| 10-01 | [Branch merges and scope-only families](2026-10-01_basidio-merges-and-scope-only.md) | decided | net 14 gained in-lineage, 2 CAAX-only PR lost, 0 label changes |
|---|---|---|---|
| 09-26 | [Anchors, full run, cap-off](2026-09-26_basidiomycota.md) | decided | 3,269 genomes; Agaricomycotina 83.0% |
| 09-26 | [Cryptococcus SXI slot](2026-09-26_cryptococcus-sxi.md) | decided | 243 genomes, 0 lost, 53 gained |
| 09-26 | [Pucciniomycotina, Wallemia, Rhodotorula](2026-09-26_pucciniomycotina-wallemia-rhodotorula.md) | decided | Wallemia 51/51; Rhodotorula P/R 61/62 held-out |
| 09-27 | [Receptor (B/PR) loci](2026-09-27_receptor-loci.md) | decided (review later) | CAAX: 6/9 vs 0/25; Agaricales PR 14->132 |

## Ascomycota (incl. Serinales, Dothideomycetes)
| Date | Report | Status | Key numbers |
|---|---|---|---|
| 09-26 | [Flank-carried rule](2026-09-26_flank-carried-rule.md) | decided | 39-bit floor: real 7/8, noise 0/52 |
| 09-26 | [Polish cap](2026-09-26_polish-cap.md) | decided | compute 135.3->74.6 h |
| 09-26 | [Serinales, C. albicans, C. auris](2026-09-26_serinales-candida.md) | decided | C. albicans 8/11 collapsed |
| 09-26 | [Dothideomycetes curation](2026-09-26_dothideomycetes.md) | decided | 90/99; C. kikuchii homothallism candidate |
| 09-27 | [Confidence rules](2026-09-27_confidence-rules.md) | decided | 86 risers; no fallback floor |

## Validation and held-out sets
| Date | Report | Status | Key numbers |
|---|---|---|---|
| 09-28 | [MAT-gene gate and validation](2026-09-28_mat-gene-gate-and-validation.md) | decided | typing 0/540 wrong; gate 96/108 vs 9/189 |
| 09-28 | [Held-out sets](2026-09-28_heldout-sets.md) | superseded by 10-01 rerun | Jena 61/64; LCG 536/621 (clean 449/533); Zygo saturated |
| 10-01 | [Held-out rerun and curator names](2026-10-01_heldout-rerun-and-curator-names.md) | decided; runtime check running | LCG 536/621, both 36->24, labels 11/1/4; Jena 61/64, both 3->1; Jena 59/64 named, LCG 61 renamed |
| 10-02 | [LCG overrides for manual review](2026-10-02_lcg-overrides-manual-review.md) | decided; manual review done 10-03 (accepted) | 12/13 rDNA agree with tree; 13 overrides (12 confirmed, 1 unconfirmed); UFBoot 75 held, 26 weak (25 inside genus only); LCG species 200->199 |
| 09-30 | [Cap V3 and regression check](2026-09-30_cap-v3-and-regression-check.md) | decided (signed off) | stress test 19/20 recovered, 0 wrong (V3); panel 166 + Zygo |

## Reviews and literature
| Date | Report | Status | Key numbers |
|---|---|---|---|
| 09-28 | [Fable review and fixes](2026-09-28_fable-review-and-fixes.md) | decided | 6 major findings; fixes at 3aec88b; 247->258 genomes |
| 09-29 | [Homothallism literature](2026-09-29_two-idiomorphs-and-homothallism.md) | decided | Z. heterogamus one locus 5.3 kb; Mycotypha ~150 kb; Syzygites two loci |
| 09-28 | [Subloci literature](2026-09-28_subloci-literature.md) | decided | one locus with subloci |

## Data quality
| Date | Report | Status | Key numbers |
|---|---|---|---|
| 09-26 | [Data quality](2026-09-26_data-quality.md) | decided | 9 amplicons suppressed |
