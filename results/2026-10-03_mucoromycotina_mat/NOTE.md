# Mucoromycotina MAT campaign (2026-10-03)

Study report: `analysis/2026-10-03_mucoromycotina-mat-campaign.md`.

## Inputs
- `inputs.tsv`: 974 genomes. BFD Mucorales + Umbelopsidales 289 (`bfd_kept.tsv`:
  288 after the suppress list), LCG 621, Jena 64. LCG and Jena are held-out sets:
  never used for training, curation or classifier builds.
- `chunk_{0..3}.tsv`, `run_chunk.slurm`: SLURM array 29371553. Code: frozen
  worktree `run-7c7ed99` (= v0.6.0). Every genome: `detect --phylum Mucoromycota
  --emit-cds-fasta`, no taxid (LCG/Jena have none; uniform mode for all).

## Detection
- `reports_all.tar.zst`: 973 `detection_report.yaml` + `wall_seconds`. The
  per-genome GFF3 and FASTA (`runs/`, 2.2 GB) are not committed.
- `collect.py` writes `calls.tsv` (one row per locus, or per uncalled genome),
  `core_proteins.faa` (894), `flank_proteins.faa`.
- 865 of 973 genomes called; 911 loci (Plus 482, Minus 425, undetermined 4).
  BFD 256/288, LCG 548/621, Jena 61/64. 47.4 CPU-h; median 183 s per genome.

## Locus size (`locus_size/`)
- Size = inner ends of the flanking genes (curator ruling 2026-10-03). Only
  polished flank models count, and only with flanks on both sides of the core
  gene: 579 loci in 23 genera. Held-out (`hold out`) and unconfirmed-override
  genomes are dropped.
- The flank pair differs by lineage (tptA|rnhA most genera; btbA|rnhA Rhizopus;
  algA|glrA Backusella, Absidia; glrA|tptA Umbelopsis). `flank_pairs` in
  `locus_size_by_genus.tsv` lists the pairs per genus.
- Medians (bp, Plus / Minus): Umbelopsis 13,063 / 12,256; Absidia 9,611 / 9,249;
  Phycomyces 6,531 / 4,166; Backusella 6,113 / 6,057; Pilaira 4,726 / 2,968;
  Rhizopus 2,619 / 1,500; Mucor 1,693 / 1,782.
- Figure `locus_size_tree.{png,svg}`: genus tree pruned from the nf_phyling
  mucoromycota_odb12 FastTree tree (reference tips preferred); Amylomyces has no
  tip.

## Synteny (`synteny/`)
- 26 loci, Plus and Minus for 13 genera (`selection.tsv`). Region = locus
  +/- 3 kb. Gene models from the genome annotation (funannotate GenBank for
  LCG/Jena; BFD GFF3); a MATPredict model is added (`*`) only where no annotated
  CDS overlaps it.
- `gbk.tar.zst` (26 GenBank files), clinker output (`clinker_alignments.txt`,
  `clinker_matrix.csv`, `mucoromycotina_MAT_synteny.clinker.html.zst`), static
  figure `mucoromycotina_MAT_synteny.{png,svg}` (`draw_static.py`, links >= 30%
  identity).

## Gene trees (`tree/`)
- `select_tree_set.py`: ingroup sexP 487, sexM 409 (hold-outs and < 60 aa
  dropped); outgroup 190 non-MAT HMG-box proteins (paralog negatives + P1).
- `run_tree.sh`: cd-hit 0.98 per idiomorph (sexP 163, sexM 173), 0.90 outgroup;
  HMG box = hmmalign to PF00505, match columns (523 seqs); full length = MAFFT
  L-INS-i + ClipKIT kpic-smart-gap (526 seqs x 526 sites). IQ-TREE 3 ModelFinder,
  UFBoot 1000, SH-aLRT 1000, seed 20261003. `run_iq_one.sh` runs one tree and
  copies the checkpoint to `ckp_<aln>/` every 30 min.
- Full length (job 29380482, 2 h 40 min): JTT+F+R8. `draw_gene_tree.py full`:
  - sexP: 163 tips; its MRCA holds 0 sexM tips and 2 outgroup tips; UFBoot 40.
  - sexM: 173 tips; its MRCA holds 0 sexP tips and 1 outgroup tip; UFBoot 93.
  - The outgroup is not one clade. The tree is rooted on one outgroup tip, and
    132 outgroup tips fall inside the sexP + sexM MRCA, so sexP and sexM are not
    shown as sisters.
  - Outgroup tips inside sexP: Dicele1 h5, GCA_016758965.1 h6; inside sexM:
    Mycafr1 h1. All are HMG genes from Mucorales genomes away from the called
    locus.
- HMG box (job 29389343, 12 h 30 min; the first run in job 29374314 hit its
  12 h limit with nothing copied back): 523 sequences x 69 sites, LG+R6.
  `draw_gene_tree.py hmg`:
  - sexP: all 162 tips in one clade with no sexM tip and 4 outgroup tips,
    UFBoot 99 (the sexP MRCA alone: 165 tips, 3 outgroup, UFBoot 93).
  - sexM: not one clade. The largest clade without sexP holds 48 of 173 sexM
    tips (UFBoot 53); the sexP clade nests among the other sexM lineages, so
    the sexM MRCA is the ingroup MRCA (435 tips, UFBoot 69).
  - 69 sites give weak deep nodes; read the sexM split as unresolved, not as
    sexP arising within sexM. The full-length tree has sexM as one clade
    (173/173, 1 outgroup tip, UFBoot 93).
- Clade test used for both trees: the largest clade that holds no tip of the
  other idiomorph (outgroup tips allowed). Full length: sexM 173/173 (UFBoot
  93), sexP 163/163 but only with 131 outgroup tips (UFBoot 29).
- RAxML-NG check on the HMG box: job 29397060 (running).
