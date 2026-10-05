# MATPredict side projects: MAT typing from reads, and MAT-like locus discovery from populations

Date: 2026-10-05
Status: draft; curator answers of 2026-10-05 recorded below. Nothing here is implemented. Separate branch and
separate code from `detect`; no change to the curated database or classifiers.

## Context
`matpredict detect` needs an assembly. Two questions it cannot answer:

1. **Typing from reads.** Many strains have reads but no assembly (or a poor
   low-coverage one). Can we call the idiomorph by aligning reads to the
   species' known MAT idiomorphs (example: Fusarium MAT1-1 vs MAT1-2)?
2. **Discovery from a population.** In lineages with no known MAT locus
   (Batrachochytrium; some Xylariales), can we find a candidate MAT-like locus
   from many strains' reads, without annotation, on the assumption that some
   HMG-box gene acts as the MAT gene?

The curator asked for both as specifications for a separate line of work
(2026-10-05).

## Data available locally (counted 2026-10-05)
| Dataset | Strains | Reads | Known MAT |
|---|---|---|---|
| ZyGoLife LCG (`/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Assembly/input`) | 881 with local reads | paired FASTQ | MATPredict assembly calls (865/973 Mucoromycotina called) |
| Batrachochytrium dendrobatidis (`Population_Genomics/B_dendrobatidis`, `Chytrid/Bd_popgen`) | about 300-330 (726 read files) | Illumina, SRA | none known |
| Rhodotorula mucilaginosa (`Population_Genomics/Rhodotorula_mucilaginosa`) | 425 rows | yes | redPR A1/A2 (MATPredict; preprint P/R assignments) |
| Rhizopus stolonifer / R. microsporus | 114 / 49 | yes | sexP/sexM (MATPredict) |
| Clavispora lusitaniae (`Population_Genomics/C_lusitaniae`) | 97 | yes | MTLa / MTLalpha |
| Aspergillus fumigatus (`Population_Genomics/Afumigatus_Global`) | 23 read files | yes | MAT1-1 / MAT1-2 |
| Xylariales | 257 BFD assemblies; no local population read set found | - | to source |

Fusarium population read sets are not local; they would come from SRA
(selection is part of the work, not done here).

## Curator decisions (J. Stajich, 2026-10-05)
- Code home: MATPredict subcommands (`matpredict reads-type`,
  `matpredict population-scan`) on a feature branch; heavier tools (minimap2,
  KMC, Bifrost) as optional pixi features.
- Hybrid / mixed strains: simulated mixes first (reads of a MAT1-1 and a
  MAT1-2 strain merged at 50:50, 80:20, 95:5) to calibrate the "both" call and
  the depth ratios; then the A. fumigatus putative hybrid strains from
  Lofgren et al. (strain list to be supplied).
- Fusarium: Fusarium oxysporum f. sp. lactucae is the main use case; the
  curator has local data (path to be supplied). Also any F. oxysporum with a
  genome, a MAT assignment and raw reads, to confirm the same assignment from
  reads.
- Batrachochytrium: all lineages first; a MAT-like partition should cut across
  lineages rather than follow them.
- Xylaria flabelliformis NC1011 is a named target (below).

## Case: Fusarium oxysporum f. sp. lactucae (Fola) -- existing read-based typing (checked 2026-10-05)
Data (curator): `/bigdata/stajichlab/nicolel/Fola/`
- Reads: `02_fastq/NCBI_SRAs_retrieved/` (SRA runs, e.g. 20a = SRR28734943;
  94 FASTQ files) and `04_AVITI_rawoutput/Martin_Nov2025/` (AVITI, 66 FASTQ
  files); assemblies in `genomes/`.
- Existing typing (N. L.), `06_Align/Mating_Types/`
  (`00_pipeline/Align/runAlign_MatingTypes_v3.sb`): reads mapped with bwa-mem2
  to the F. oxysporum MAT1-1 and MAT1-2 idiomorph sequences (GenBank
  AB011379.2, 5,220 bp; AB011378.1, 5,122 bp; `References/`); `samtools
  coverage` per idiomorph; breadth separates the types (strain 20a: MAT1-2
  99.6% breadth, 74x; MAT1-1 11.6%, the part shared with the flanks).
- Result table `MAT_vs_Phenotype_counts.tsv`: 148 strains -- MAT2 129, MAT1 17,
  Both 1, Neither 1, cross-tabulated with race (MAT1 strains are all
  weak/non-pathogenic in that table). VSP-0947 has a de novo assembled MAT locus
  (`VSP-0947_denovo/`).
- Use for Tool A: this is the alignment path, already run; Tool A should
  reproduce these 148 calls from the same reads (the benchmark), and add the
  k-mer path, the depth ratios for "Both", and a standard report.
- MATPredict database: Fusarium records exist for F. fujikuroi (5127_mo44
  MAT1-1, 5127_mo45 MAT1-2) and F. graminearum (5518_3639 combined), none for
  F. oxysporum. Curating AB011379.2 and AB011378.1 as records would give the
  same-species reference panel and let `detect` type the Fola assemblies, so
  assembly calls, read calls and N. L.'s calls can be compared three ways.

## Case: Aspergillus fumigatus putative hybrids (Lofgren et al. 2022)
Source: Lofgren LA et al. 2022, PLoS Biol, doi:10.1371/journal.pbio.3001890
(PMC9714929); S9 Fig; supplementary table (journal.pbio.3001890.s009).
- Their method, summarised: MAT type from BLASTN of the MAT1-1 (AY898661.1,
  strain AF250) and MAT1-2 (Afu3g06170, Af293) CDS against each assembly; the
  11 strains with hits to both were rechecked by mapping raw reads to the two
  references (Bowtie2, very-sensitive-local), with depth profiles, consensus
  sequences, and ploidy checks (k-mer spectra with GenomeScope, k = 21; allele
  frequencies at heterozygous sites). S9 Fig shows 9 strains with alignments
  over both idiomorphs at different depths.
- Design point for Tool A: about 270 bp at the end of the MAT1-2 reference is
  shared with MAT1-1, so a MAT1-1 strain covers that part of the MAT1-2
  reference. Breadth and depth must be computed on idiomorph-specific
  positions only (mask the shared part), as the F. oxysporum data also show
  (MAT1-1 11.6% breadth in a MAT1-2 strain).
- Putative hybrids named by the curator (reasonable read coverage on both
  idiomorphs): IFM_59359, AF100-1_3, IFM_61407.
- Local data (checked 2026-10-05):
  - assemblies: `/bigdata/stajichlab/shared/projects/Afumigatus_pangenome/scaffolded/genomes/`
    (610 files; IFM_59359 and IFM_61407 present; AF100-1_3 not found there);
  - reads: full alignments to Af293 (FungiDB-50) as CRAM in
    `/bigdata/stajichlab/shared/projects/Population_Genomics/Afumigatus_Global/aln/`
    (IFM_59359 0.53 GB, AF100-1_3 1.49 GB, IFM_61407 0.49 GB); use these
    (samtools fastq, or the alignments directly);
  - do NOT use `Afumigatus_Global/unmapped/*.fastq.gz`: they are only the reads
    that did not map to Af293 (input to `pipeline/05_assemble_unmapped.sh`), so
    the Af293-type (MAT1-2) reads are missing from them.
- Test: Tool A must report "both" with depth ratios for these 3, single
  idiomorphs for the population, and match the published calls; the ploidy
  check (k-mer spectrum, allele balance) is reported next to the call, not used
  to make it.

## Results 2026-10-05: what `detect` (assembly) says for the two read-typing test sets
- A. fumigatus putative hybrids (`results/2026-10-05_afum_hybrids_detect/`,
  main at the time; assemblies `/bigdata/stajichlab/shared/projects/Afumigatus_pangenome/genomes/`,
  symlinks into `Afum_popgenome/asm/scaffold/genomes_scaffolded/`):
  - BLASTN of MAT1-1 (AY898661.1) and MAT1-2 (Afu3g06170 from Af293 FungiDB-50)
    finds both in all three: MAT1-1 at 99.7-99.8% over 2.3-2.4 kb on a small
    scaffold (IFM_59359 scaffold_79, 2,420 bp; DMC_AF100-1_3 scaffold_134,
    2,318 bp; IFM_61407 scaffold_69, 2,418 bp); MAT1-2 at 99.1-99.4% over
    0.9-1.1 kb on scaffold_3 (the chromosome-3 scaffold).
  - `detect` calls only MAT1-2 (medium, mat_locus on scaffold_3 with SLA2,
    APN2, COX13, MAT1-2-4). The MAT1-1 contig is withheld as
    `below_fraction_floor` (DMC_AF100-1_3, IFM_61407) or not reported
    (IFM_59359): a 2.4-kb contig carries no flank genes.
  - So the assembly path under-reports these hybrids; read depth (Tool A) is
    what separates a second nucleus/allele from contamination. This is the
    case for the report-only "unlinked second idiomorph" check.
  - DMC_AF100-1_3 and DMC2_AF100-1_3 are two near-identical scaffoldings (872
    scaffolds each, sizes differ by 180 bp); BLAST results are the same.
- Fola assemblies (`results/2026-10-05_fola_detect/`, run-94c1a3b with the new
  F. oxysporum records): 4 of 19 typed -- AT141 MAT1-2 (Chr7), JCP043 MAT1-2
  (chr8), VSP-0916 flye MAT1-2, VSP-0980 MAT1-1 (all high, mat_locus). The
  other 15 `*.fna` files in `/bigdata/stajichlab/nicolel/Fola/genomes/` are not
  readable by the jstajich account (Permission denied); they need group read
  permission before a re-run. Tool bug found: an unreadable genome surfaces as
  a BLAST "No alias or index file" error instead of a clear input error.

## Case: Xylaria flabelliformis NC1011 (checked 2026-10-05)
What MATPredict v0.6.0 found (`results/2026-10-03_ascomycota_v060/`, wave_7):
- GCA_022453505.1 (JGI Xylcub1, NC1011): one call, Ascomycota:MAT MAT1-2,
  medium, idiomorph_gene_only, JAJLYR010000012.1:276,138-285,941. The gene is
  KAI0192626.1 ("HMG box protein", 734 aa; PF00505 E 5e-23), matched only
  weakly by MAT1-2-1 (44%) and MAT1-1-3 (41%) references, among housekeeping
  genes (dynactin, SRP19, TIP49, TFIID). Its length and neighbours suggest a
  general HMG transcription factor, not MAT1-2-1 (not tested further).
- The conserved flank block is intact on JAJLYR010000004.1: APC5
  (KAI0195599.1, 240 kb) - SLA2 (KAI0195601.1, 244,382-247,733) - COX13
  (KAI0195602.1) - APN2 (KAI0195603.1, 249,891-252,058), with no MAT gene
  between SLA2 and APN2. MATPredict withheld this cluster (no MAT gene).
- Strain G536 (GCA_007182795.1) shows the same pattern (flank block on
  VFLP01000025.1, no call).
- So no MAT locus is identified: the idiomorph genes are either moved away
  from SLA2-APN2 (as in A. nidulans), too divergent for the current
  references, or absent. Next steps: search the proteome for alpha-box
  (PF04769) and HMG proteins genome-wide and rank by similarity to
  Sordariomycetes MAT proteins; Tool B if a Xylaria population read set exists.
- Further work on NC1011 (interval, gene-order classes, HMG tree, RNA-seq):
  branch `xylariales-interval`, `docs/HANDOFF-xylariales-2026-10-05.md`.
- SRA (checked 2026-10-05, txid2512241): no population. DNA from two strains
  only: NC1011 PacBio Sequel WGS (SRR8861568-73, ~12 Gb; the JGI assembly) and
  G536 Illumina MiSeq WGS (SRR9166620, 14 Gb, filed as "Xylaria cubensis").
  NC1011 RNA-seq: transcriptome SRR8861595 (21 Gb) and 7 expression-profiling
  runs (SRR37043446-52, 2-3 Gb each). Use: expression evidence and gene-model
  checks for candidate alpha-box/HMG genes in NC1011; Tool B is not possible.

## Tool A: idiomorph typing from reads (`matpredict reads-type`, proposed)

### Input and reference
- Reads: paired or single FASTQ, or an SRA run accession; a taxid or species.
- Reference panel, built from the curated database and MATPredict calls in
  the same species (or genus):
  - idiomorph-specific sequence of each idiomorph (core genes and the region
    between the inner flank ends);
  - shared flank genes on both sides (controls: they must be covered in every
    strain);
  - optional: a single-copy genome control set (for example BUSCO genes) for
    depth normalisation.

### Method
1. Fast path (k-mers): k-mers unique to each idiomorph (absent from the other
   idiomorph and from the flanks); count their presence in the reads
   (streaming counter or KMC). Gives a call in minutes and works at very low
   depth.
2. Alignment path: map reads (minimap2 `-x sr` or bwa-mem2) to the panel;
   per region, breadth of coverage and depth normalised to the flanks.
3. Divergent references (genus level, no same-species record): translated
   search of reads against the idiomorph proteins (DIAMOND blastx), reported
   as lower confidence.
4. Call per strain: one idiomorph; both (homothallic, diploid/heterokaryon,
   or mixed sample — reported with depth ratios, not resolved); none (reference
   gap or too divergent); with breadth, depth ratio and k-mer counts as evidence.

### Validation
- LCG: 881 strains with reads and assembly-based calls. Concordance of the read
  call with the assembly call, by genus; discordant cases reviewed (sample
  mixup and misnaming are known in this set).
- Downsampling: 0.5x, 1x, 2x, 5x, 10x genome depth on 50 strains: the lowest
  depth that keeps concordance >= 95% (to be measured).
- Controls with known MAT: A. fumigatus, C. lusitaniae, R. mucilaginosa,
  Rhizopus; Fusarium oxysporum f. sp. lactucae (curator's local data) and other
  F. oxysporum with a genome, MAT assignment and reads; homothallic
  F. graminearum should give "both".
- Mixed and hybrid samples: simulated mixes of a MAT1-1 and a MAT1-2 read set
  (50:50, 80:20, 95:5) first, then the Lofgren et al. A. fumigatus putative
  hybrids.
- Zygo 23 truth set where reads exist.

### Outputs
Per strain: call, confidence, evidence (breadth, depth ratios, unique-k-mer
counts per idiomorph), reference panel version. Per run: concordance table and
a summary for the aggregate report.

## Tool B: MAT-like locus discovery from a population (proposed)

### Idea
In a heterothallic species, the MAT locus is a balanced polymorphism: two (or
more) mutually exclusive haplotype blocks (idiomorphs), each present in part of
the population, between flanks present in all strains. Its pattern is unusual:
the partition of strains is often not explained by the strain phylogeny, and
the blocks carry a transcription-factor gene (HMG box, alpha box or
homeodomain). Search the population for regions with that pattern.

### Method
1. Reference-based presence/absence:
   - map all strains to one reference assembly (and, when known, a second of
     the other type); coverage in windows (for example 1 kb), normalised per
     strain;
   - keep windows whose normalised depth is near 0 in some strains and near
     the genome level in others (bimodal), in 10-90% of strains;
   - group adjacent windows with the same strain pattern into blocks;
   - for strains lacking a block, assemble their reads that do not map to the
     reference near the block's flanks, to recover the complementary
     (other-idiomorph) sequence.
2. Reference-free k-mers (alternative and check):
   - per-strain k-mer sets; unitig presence/absence matrix (for example
     KMC + Bifrost);
   - find unitig sets with mutually exclusive presence (set A in group G, set B
     in its complement), assemble each set, annotate.
3. Ranking each candidate block:
   - partition balance (MAT is often near 50:50; skewed pops are expected);
   - independence from the strain phylogeny (consistency index / homoplasy of
     the partition on a core-genome tree; MAT partitions recur across clades,
     accessory regions usually do not);
   - gene content: HMG box (PF00505), alpha box, homeodomain, pheromone /
     receptor genes (hmmscan on translated contigs);
   - flank genes present in all strains on both sides.
4. Diploids (Bd is diploid): idiomorph states are 0, 0.5 and 1 relative depth
   (absent, heterozygous, homozygous); loss of heterozygosity in Bd lineages
   adds 0/1 patterns that are not MAT. The model must use the three-state depth.

### Evidence levels for a candidate (proposed; thresholds calibrated on the controls)
Presence/absence has many non-MAT causes (transposons, accessory chromosomes,
deletions, contamination, assembly gaps), so a candidate needs several
independent signals:
1. Clean presence/absence: near-zero depth over the whole block in some
   strains, present in others; minority group >= 3 strains and >= 10%; the
   flanking sequence on both sides present in every strain.
2. Complementary block: strains lacking block A carry a different sequence B
   between the same flanks (assembled from their unplaced reads or k-mers):
   the idiomorph signature.
3. Gene content: A or B carries a transcription-factor gene (HMG box, alpha
   box, homeodomain).
4. Independence from the strain tree: the A/B partition recurs in several
   clades rather than marking one lineage.
5. Artefact checks: not repeat-dominated, not at a contig or chromosome end,
   not a whole-chromosome depth change.
Levels: strong = 1-4 with 5 passing; weak = 1 + 3, or 1 + 2; otherwise
reported only as a presence/absence polymorphism. Each control's known MAT
locus must reach "strong"; the number of non-MAT strong candidates per control
is the false-positive measure.

### Positive controls (run blind; the known locus must rank near the top)
A. fumigatus (MAT1-1/MAT1-2), C. lusitaniae (MTL), R. mucilaginosa (redPR/redHD),
Rhizopus stolonifer and R. microsporus (sexP/sexM). Report the rank of the known
locus and the false candidates above it.

### Targets
- Batrachochytrium dendrobatidis: about 300-330 strains locally; no MAT locus
  known; diploid; lineages largely clonal (risk: partitions that follow
  lineages).
- Xylariales: 257 BFD assemblies; a species with enough strains (>= 30) and
  reads must be found first.

### Outputs
Candidate table (block coordinates, strain partition, balance, phylogenetic
independence, genes and domains, complementary contigs), figures (presence
matrix ordered by the strain tree; block gene order), and a short report per
species.

## Risks and limits
- Presence/absence also comes from transposons, accessory chromosomes,
  aneuploidy, contamination and uneven coverage; ranking must not rely on one
  signal.
- Clonal populations make a MAT partition look like a lineage marker.
- Homothallic or asexual species may have no partition at all.
- Read typing depends on the reference panel; a strain with an idiomorph
  absent from the panel will read as "none".

## Plan (phases)
| Phase | Content | Done when |
|---|---|---|
| A1 | Reference-panel builder from db records + same-species calls | panels for 10 species |
| A2 | `reads-type` k-mer and alignment paths; per-strain report | runs on LCG subset |
| A3 | LCG concordance + downsampling; controls | concordance table, depth floor |
| B1 | Window coverage + bimodality + blocks on a control species (C. lusitaniae, haploid) | known MTL ranks first or near |
| B2 | Phylogenetic-independence score; k-mer path; other controls | known loci recovered across 4 controls |
| B3 | Bd run; Xylariales data sourcing and run | candidate reports |

Compute (not measured): mapping ~300 strains to one reference is in the range
of a few hundred CPU-hours; the k-mer path needs per-strain k-mer databases on
`$SCRATCH`.

## Questions for the reviewer
Answered 2026-10-05: code home, hybrid order, Fusarium source, Bd scope (see
Curator decisions). Still open:
1. Accept the proposed evidence levels for Tool B candidates (above)?
2. Which Xylariales species has a population read set (>= 30 strains)? (X. flabelliformis: none in SRA.)
3. A. fumigatus AF100-1_3 assembly location (reads are in the CRAM set).
4. F. oxysporum MAT1-1/MAT1-2 records: curator said yes (2026-10-05); curated
   on branch curate-foxysporum.
