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
1. Minimum evidence to report a "candidate MAT-like locus" from Tool B.
2. Which Xylariales species has a population read set (>= 30 strains)?
3. Paths: F. oxysporum f. sp. lactucae data; the Lofgren et al. hybrid strain
   list.
