# Flank-first MAT scan: architectures in lineages with no defined MAT
Status: stage 1 done on 496 BFD genomes; stage 2 (RBH flank collinearity) run on Rhizophagus only and did not recover the control

## Question
Do lineages with no curated MAT locus carry gene architectures like known MAT loci?

## Method (stage 1)
- Genes from funannotate GFF3 and proteins (BFD genome_annotation). Window = 5 genes each side.
- Candidates: Pfam 38.2 GA-cutoff hits to HMG_box, HMG_box_2, Homeodomain, and diamond hits (E<=1e-5) to controls.
- Controls: 93 Ropars 2016 proteins (HD1-like, HD2, choline transporter-like, phosphoglycerate mutase; PMID 27572831) and the 41 curated Mucoromycota MAT proteins (sexM, sexP, tptA, rnhA, glrA, algA, btbA).
- Chance baseline: expected number of genomes with an adjacent pair if the genes sat at random (see scripts/flank_first_null.py). It ignores contig breaks and so over-states chance adjacency.
- Code: scripts/flank_first_scan.py, flank_first_summarize.py, flank_first_null.py, flank_first_pairs.py (branch flank-first-scan). Jobs 32062718 (stage 1).
- Results: /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-13_flank_first_scan/stage1_all/ (summary_by_lineage.tsv, null_by_lineage.tsv, pairs_by_lineage.tsv, per_genome/).

## Results
1. Control (Rhizophagus): HD1+HD2 within 5 genes in 27 of 29 genomes. In 31 of 31 pairs the two genes are divergent (-/+) and 0 genes apart. Chance expectation is 0.0.
2. Other Glomeromycota (30 genomes, 14 genera): 9 have the same adjacent divergent pair; 8 after the contamination check in 2026-11-13_amf-hd-locus.md, which removes Oehlia diaphana SD1 (its locus contig is a copy of Rhizophagus sp. MUCL 43196 sequence) (Ambispora 2/2, Diversispora 2/2, Glomus 2/2, Entrophospora 1/3, Geosiphon 1/1, Oehlia 1/1). In 12 more, HD2 or HD1 is present but not as a pair (Acaulospora and Gigaspora: both genes, on different contigs; Cetraspora, Dentiscutata, Paraglomus, Racocetra, Scutellospora, Entrophospora: HD2 only). Funneliformis (0/4) and one Paraglomus have no hit. Hits are strong (E 1e-30 to 0).
3. Mortierellomycota (100 genomes): no HD1/HD2 hit, no adjacent Homeodomain pair, no HMG+Homeodomain pair, no sexM/sexP+flank locus. Chance expected 3.9, 11.4 and 3.3 genomes respectively.
4. Kickxellomycota: HMG+Homeodomain adjacent in 71 of 177 genomes (27.1 expected by chance); 63 of the 71 are Coemansia.
5. Basidiobolomycota: two adjacent Homeodomain genes in 15 of 18 genomes (5.3 expected), all Basidiobolus. All 19 pairs are convergent, not divergent.
6. Chytridiomycota: isolated divergent adjacent Homeodomain pairs in 10 of 83 genomes (Chytriomyces 9 pairs, Entophlyctis 4); several genomes have two pairs. Other chytrid homeodomain adjacency comes from large same-strand arrays (1,207 pairs in 28 genomes). sexP/sexM+rnhA/glrA within 5 genes in 4 genomes (3 Chytriomyces hyalinus, 1 Boothiomyces), against 2.0 expected.
7. Zoopagomycota (11 genomes): 3 genomes (Cochlonema, Piptocephalis, Thamnocephalis) have two adjacent Homeodomain genes; 2 of 3 are convergent. Blastocladiomycota (14): nothing above chance.

8. Known flank genes next to the HD pair (stage 1, window 10 genes). In Rhizophagus genomes with an HD hit (28), a choline transporter-like gene is within 10 genes in 26 and a phosphoglycerate mutase gene in 23. Chance expectation is 0.05 and 0.02 genomes. In other AMF with an HD hit (22), choline 11 (0.04 expected) and PGAM 4 (0.01 expected). Genera with the choline neighbour: Ambispora 2/2, Diversispora 2/2, Glomus 2/2 (plus PGAM 2/2), Geosiphon 1/1 (plus PGAM), Oehlia 1/1, Paraglomus 2/2 (PGAM 1), Gigaspora 1/3. None in Acaulospora, Cetraspora, Dentiscutata, Entrophospora, Racocetra, Scutellospora. The chance model is the single-copy approximation 1-(1-21/N)^n; the choline-like control is a gene family (all 100 Mortierellomycota genomes have a hit), so these odds are approximate.
9. Stage 2 (job 32062611, Rhizophagus, reference DAOM-197198, RBH E<=1e-10, 5-gene window). It did NOT recover the control locus. The 31 HD1+HD2 loci have collinear blocks of 0 to 5 flank genes in 28 of 31 (3 loci reach 10); the null percentile median is 0.21. Across all loci, 42% have a block >= 6. The metric therefore does not separate the control locus and I do not use it to rank candidates. Likely causes (not tested): fragmented assemblies, 31,000 to 52,000 gene models per genome, and one reference for 29 isolates. Flank-first search with known flank genes (item 8) worked; unbiased flank conservation did not.

## Limits
- Hits to HD1/HD2 show sequence similarity to the AMF proteins. They do not show a mating function. The authors call the Rhizophagus locus putative; no source I read tests it.
- sexM/sexP diamond hits occur in almost every genome (HMG-box family). They are not specific.
- The Basidiobolus, Coemansia and Chytriomyces patterns are genus-wide. They could be conserved gene pairs unrelated to mating.
- No phenotype or mating data were used.
- Stage 2 (RBH flank collinearity) did not recover the control locus; see item 9.
