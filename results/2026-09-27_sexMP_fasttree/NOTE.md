# sexM/sexP tree, FastTree redraw on the classifier-scan labels (2026-09-27)

- Loci: Mucoromycota from results/2026-09-26_hmm_classifier/Mucoromycota_1a00b0a
  (623 sexM/sexP loci; 35 genomes changed, 26 new locus ids, re-extracted);
  Mortierellomycota/Kickxellomycota carried over (1,145 loci). 1,768 loci total.
- Same outgroups (non-locus HMG copies minus any a new locus covers; 2 MAT1-2-1),
  same hmmalign PF00505 match columns (69), same dedup: 868 tips.
- FastTree 2.1.11 `-lg -gamma`. Support = FastTree local support (0-1), not UFBoot;
  no orthology is claimed from it. sexP clade: all 6 refs, local support 0.23
  (IQ-TREE UFBoot was 99). sexM: 7 of 9 refs in the best pure group (0.76).
- Rooted on Fusarium graminearum MAT1-2-1 (5518_3639_MAT_combined).
- Concordance, called loci placed in sexP or sexM (concordance.txt):
  before 148/195, after 204/206. Remaining called disagreements
  (disagreements.tsv): Umbelopsis sp. M5902 (Minus, low, no classifier call;
  HMG box in sexP clade) and Mucor hiemalis gzMucHiem1 (Plus, low, no classifier
  call; HMG box in sexM group). One withheld: Endogone sp. FLAS-F59071.
- Files: tree_ft.pdf, tree_ft.treefile, tree_ft.named.{treefile,nex},
  tip_names.tsv, candidates.tsv, concordance.txt, disagreements.tsv.
