# HMM idiomorph classifier (Mucoromycota:MAT) and the all-genes cap rank

Code: polish-scope-cuts dfa31fa (cap rank), 7277e10 (classifier), 67a122c,
6e2f58d, 1a00b0a (fixes found by the 293-genome runs). Final runs on 1a00b0a.

## Classifier build (db/Mucoromycota/classifiers/MAT, scripts/build_idiomorph_hmms.py)
Training (curator option (a)): curated sexM (9) / sexP (6) records + 70 UFBoot-99
sexP-clade proteins (make_training_extra.py; disputed sets and Zygo excluded).
sexP 76 seqs / 75 unique / 27 genera, identity 14.1/30.1/100; sexM 9/9/7 genera,
19.0/29.7/50.0. Leave-one-genus-out 85/85, worst correct margin 22.5 bits.
min_margin 11.2 = half that.

## Held-out performance
- Zygo 23 proteins scored directly (score_heldout.py): 23/23, smallest margin 47.5.
- Zygo 23 pipeline, 1a00b0a: scaffold 23/23 locus, 23/23 idiomorph (classifier
  decided 11, fallback 12 -- annotated proteome fast path, nothing polished);
  contig 23/23, 23/23 (classifier 22, fallback 1). Smallest classifier margin 34.0.
- 49 disputed proteins (heldout_scores.tsv): 13/13 Minus-in-sexP-clade -> sexP
  (margins 50.8-236.9); 33/33 Rhizopus sexM-only -> sexM (-109.5 to -169.7);
  Lichtheimiaceae/Syncephalastraceae 2 sexP (135.7, 144.4), 1 at -2.3 (below
  min_margin: undetermined).
- 2 Umbelopsis Minus calls: U. ramanniana AG -> Plus (81.6 vs 30.6, margin 51.0).
  Umbelopsis sp. M5902 stays Minus on the fallback: no core protein was modelled.

## 293-genome re-scan vs 4a24ffb (compare_output.txt, label_changes.tsv)
238 unchanged, 4 Plus->Minus, 1 Minus->Plus, 2 Minus->undetermined, 0 lost, 0 new.
- Plus->Minus: Mucor irregularis B50 and B7584 (Minus 80.7 vs Plus 50.6),
  Mucor ardhlaengiktus CBS 210.80 (77.2 vs 53.0), Syncephalastrum racemosum
  NRRL 2496 (64.4 vs 52.4; margin 12.0, just above min_margin).
- Minus->Plus: Umbelopsis ramanniana AG.
- undetermined: Syncephalastrum monosporum B8922 and PYS2302 (65.2 vs 54.9).
Calls decided by the classifier 232/245 (decisive margins min 11.7, median
169.7); 13 fall back (no modelled core protein).
Runtime: median wall per genome 119 -> 128 s; total 9.65 -> 9.86 h.

## Cap rank (compare_cap.py, cap_compare.txt; 34 genomes, Dothideomycetes records)
All-genes rank vs live rank: 0 lost, 2 gained, 0 changed. Gained: Acidiella
bohemica MAT1-2 high, Bauco1 MAT1-2 medium. Zymoseptoria brevis was called by
both runs. Aulographum stays uncalled: its best cluster (83.6%) has 2 genes and
still ranks outside the top 6. GCA_029290875.1 (the predicted loss) had no MAT
call in either run. Median wall 258 vs 270 s.

## Bugs found by the real runs and fixed
- 7277e10 flipped every loser hit in a cluster: 11 real calls lost (67a122c).
- Models off frame by one base translated to junk: 3 Mucor circinelloides
  undetermined at ~0 bits (6e2f58d, 3-frame translation).
- Pooled verdict applied to a second HMG gene: 2 calls lost (1a00b0a, per-position).

## Limits
Label changes have no mating-type truth outside Zygo. The sexM training set is
9 records from 7 genera; sexM clade is not supported in the tree.
