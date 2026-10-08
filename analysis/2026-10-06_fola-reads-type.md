# Fola: MAT typing from reads only (`matpredict reads-type`), plus assembly checks
Status: open (first version built and benchmarked; no independent test yet)

## Question
Can the MAT idiomorph of a Fusarium oxysporum f. sp. lactucae (Fola) strain be called from its raw reads alone,
and does the call match the existing samtools-breadth calls (N. L.) and the assembly calls (`detect`)?

## Data and code version
- Code: branch `reads-type` (worktree `.claude/worktrees/reads-type`), new package `src/MATPredict/reads/`,
  18 tests in `tests/reads/`. `detect` runs use frozen worktree `run-94c1a3b` (adds the F. oxysporum records AB011379.2 / AB011378.1).
- Panel: AB011379.2 (MAT1-1, 5,220 bp) and AB011378.1 (MAT1-2, 5,122 bp) from N. L.'s `06_Align/Mating_Types/References/`.
- Reads: 148 strains, trimmed pairs in `/rhome/nicolel/bigdata/Fola/05_Pre-Processing/`.
- Truth: N. L.'s `06_Align/Mating_Types/Coverage/MAT_coverage.tsv` (bwa-mem2 per-idiomorph breadth, now recorded in the CRAM headers) with the rule breadth >= 90%.

## Method
1. Step 1, assembly path: `detect --taxid 5507` on 19 Fola assemblies plus the VSP-0947 SPAdes scaffolds (job 29546488).
2. Step 2, isolate 50a: whole-genome SPAdes `--isolate` (job 29546493), BLASTN of both references, `detect` (job 29548342).
3. Step 3, `reads-type`: k = 31. Unique k-mers of each idiomorph are those absent from the other idiomorph. K-mers present in both (the shared flank
   ends, 282 k-mers) are a single-copy depth control. Reads are matched against both strands. First 8,000,000 reads per strain (R1 first). Per idiomorph: breadth
   (fraction of unique k-mers seen) and depth (mean count).
   Present = depth >= 0.10 x shared depth AND breadth >= 0.50 x (1 - exp(-depth)). Call: one idiomorph, `both`, `none`, or `low_depth` (shared depth < 1).
   `trace_<idiomorph>` flags a non-zero signal that fails the rule.

## Results
Step 1 (`results/2026-10-06_fola_detect_all/runs/*/detection_report.yaml`):
- All 19 assemblies are typed with one high-confidence `mat_locus`: 17 MAT1-2 and 2 MAT1-1 (FON_AJ275, VSP-0980). Fo47 and Fo47_REF_assembly are the same genome (two files, both MAT1-2); counted by file, not by genome.
- VSP-0947 SPAdes scaffolds: BLASTN finds both idiomorphs, each on a short contig. MAT1-1: 5,220 bp at 99.9% on NODE_931 (5.5 kb, cov 9.7). MAT1-2: 4,073 bp at 99.4% on NODE_1613 (4.1 kb, cov 7.7). `detect` calls only MAT1-1 (medium, `idiomorph_gene_only`), because neither contig has flanking genes.
- Every assembly also has 6 to 10 withheld `modelled_gene_bar` loci. These are the same low-quality background that all genomes show; none was examined.

Step 2, isolate 50a (`results/2026-10-06_fola_50a/`):
- MAT1-2 locus complete on NODE_86 (169.7 kb, cov 30.7): 99.1% over 5,106 bp. `detect`: MAT1-2 high, `mat_locus`, NODE_86:116,422-128,495, genes COX13, APN2, MAT1-2-1, SLA2.
- MAT1-1 shows only as fragments of 362 to 460 bp at 86 to 90% identity, on contigs of coverage about 1 (against 30.7). There is no MAT1-1 idiomorph in the assembly.
- samtools gave 44.5% MAT-1 breadth for 50a. Only 20 reads lie in the middle of the MAT-1 reference, with 5 to 17 mismatches per 150 bp. The k-mer path gives MAT1-1 depth 0.25, the same as clean MAT-2 strain 20a (0.25).

Step 3, concordance with the samtools call (`results/2026-10-06_fola_reads_type/v2_concordance.tsv`, n = 148):
| Truth | k-mer call | n |
|---|---|---|
| MAT1-1 | MAT1-1 | 17 |
| MAT1-2 | MAT1-2 | 127 |
| MAT1-2 | low_depth | 2 (AL185, VSP-0992) |
| both | both | 1 (VSP-0947) |
| none | low_depth | 1 (VSP-0931) |
- Agreement 145/148 (98.0%). The 3 differences are refusals, not wrong calls. No strain got the opposite idiomorph.
- AL185 and VSP-0992 have shared depth 0.78 and 0.83, just below the floor of 1, while MAT1-2 depth is 2.64 and 0.63.
- VSP-0931 has no signal in either method (samtools MAT-2 breadth 0.98%, k-mer breadth 0).
- Speed: about 1.5 min per strain at 8M reads in pure Python on 1 CPU.

## Limits
- The 50% breadth fraction was set after a first run with a 90% fraction failed. MAT1-2 breadth was 0.78 of expected in all strains, because the Fola MAT1-2 allele differs from the GenBank reference at about 22% of the unique k-mers. The new value was chosen from the failed first run (the first 40 strains it finished, all MAT1-2 at 0.78 of expected) and a k scan on 20a and VSP-0947. The same 148 strains then gave the 98.0%, so it is a consistency check on the same cohort, not an independent test.
- A test with a SNP every 40 bp (about half the k-mers lost) is called `none`; every 60 bp is called correctly. Allele divergence above about 1 SNP per 40 bp needs the alignment path or a shorter k (not built).
- The truth comes from the same reads, so both methods share any sample mix-up. Only one "both" strain (VSP-0947) tests the non-single outcome, and the one truth "none" strain (VSP-0931) came out `low_depth`, so `none` is not tested on real data (only in a synthetic test).
- Mixed samples are not calibrated yet. The planned simulated mixes (50:50, 80:20, 95:5) are not run. A 2% contaminant (like 50a's trace) would fall under the 0.10 depth cut-off and be flagged, not called. (Superseded in part by Update 4: mixes were run for the k-mer path, 3 strain pairs.)
- Only the first 8M reads (R1 first) were used. Depth is a k-mer count, not a genome depth.
- Alignment path, divergent references (blastx) and a panel builder from the database are not built.

## Update: panel built from Fola genomes (same day)
Why: the first panel (GenBank) lost about 22% of the MAT1-2 unique k-mers in every strain. The Fola locus is 99.119% identical to AB011378.1 over 5,106 bp
(39 mismatches, 6 gaps), with the same alignment in AT141, JCP043 and the VSP-0916 flye assembly. MAT1-1 is 99.789% (11 mismatches) in VSP-0980.
Panel (`panels/fusarium_oxysporum_fola/`, README there): MAT1-2 from AT141 Chr7:1107520-1112622 (chromosome-level Fola genome), MAT1-1 from VSP-0980 NODE_12:304601-309820.
FON_AJ275 was not used: it is f. sp. niveum, not Fola.
Run: `run_v3.slurm`, same code, same 8M-read cap, same rule (job 29550029). Result: `results/2026-10-06_fola_reads_type/v3_concordance.tsv`.

| Panel | agree | MAT1-2 strains: median MAT1-2 breadth (min) | MAT1-1 strains: median MAT1-1 breadth (min) |
|---|---|---|---|
| GenBank (v2) | 145/148 (98.0%) | 0.783 (0.517) | 0.935 (0.704) |
| Fola (v3) | 147/148 (99.3%) | 1.000 (0.627) | 0.975 (0.732) |

- v3 table: MAT1-1 17/17, MAT1-2 129/129, both 1/1. AL185 and VSP-0992 (low_depth in v2) are now MAT1-2. The one difference is VSP-0931 (truth none, `low_depth`): no signal in either method.
- 50a: MAT1-2 breadth 1.000, MAT1-1 breadth 0.011 (depth 0.02), flag `trace_MAT1-1`. 20a has no MAT-1 signal at all. VSP-0947: both, breadth 0.968 and 1.000, depths 7.5 and 7.3.
- With breadth now near 1.0 where the allele matches, the 0.50 breadth fraction is no longer needed to pass the cohort; it is kept unchanged. It is untested on the Fola panel as a tighter value, so this note does not claim a better threshold.
- Circularity: VSP-0980 supplied the MAT1-1 reference, so its call (breadth 0.732) is not independent. AT141 has no reads in the 148, so the MAT1-2 reference strain is not in the benchmark. The other 147 strains are independent of the reference intervals, but the truth (samtools breadth) is still from the same reads.
- The v2 limits above stay true for a panel taken from another species or lineage.

## Update 2: no species-specific locus: protein (DIAMOND blastx) search of reads
Question (curator): without a locus defined for the species or population, can blastx of reads against MAT proteins match the k-mer accuracy?
Test design: references are database proteins (`scripts/build_blastx_tiers.py`), each tier holds one taxonomic distance from *F. oxysporum*:
T1 same species (2 proteins: MAT1-1-1, MAT1-2-1, other formae speciales), T2 same genus, other species (F. fujikuroi, F. graminearum; 8 proteins),
T3 other Hypocreales (0 proteins in the database, not tested), T4 Pezizomycotina outside Hypocreales (49 MAT proteins in 25 MAT1-1 and 24 MAT1-2 entries).
Every tier also holds 37 APN2/SLA2 proteins as a depth control. DIAMOND 2.2.6 blastx, `-k 3 -e 1e-5`, first 4,000,000 R1 reads of each of the 148 strains (job 29551121).
Call: best hit per read (max bitscore); keep reads with identity >= I and aligned length >= 30 aa; density of a class = reads of its best gene / mean reference length;
present = density >= 0.10 x control density (mean of APN2 and SLA2) and >= 3 reads. Truth = the samtools breadth call, as above (`v3_concordance.tsv`).
The filter was tuned on odd-numbered strains (sorted by name) and tested on even ones (`scripts/blastx_threshold_sweep.py`, `sweep_calls.tsv`).

| Tier | min identity | tune agree (74) | test agree (74) |
|---|---|---|---|
| T1 same species | 0 / 40-80 | 87.8% / 100% | 90.5% / 97.3% |
| T2 same genus | 0 / 40-80 | 95.9% / 100% | 97.3% / 98.6% (100% at 80) |
| T4 distant | 0 | 83.8% | 89.2% |
| T4 distant | 40 / 50 | 79.7% / 73.0% | 89.2% / 87.8% |
| T4 distant | 60 or more | 8% or less | 9.5% or less |

- Close references (T1, T2) reach 97 to 100% from a filter of 40% identity upward. The tune split cannot choose between 40 and 80%; every value gives 100%.
- At T1 (identity 50) and T2 (identity 60) over all 148 strains: MAT1-1 17/17, MAT1-2 128/129 and 129/129, both 1/1. Differences: VSP-0931 (below) and, at T1 only, VSP-0992 (MAT1-2 called none; shared control depth 0.85).
- T4 (distant references) is asymmetric. At identity 0: MAT1-2 strains 127/129 correct; MAT1-1 strains 1/17 correct (9 called MAT1-2, 2 called both, 5 none). The false MAT1-2 calls in MAT1-1 strains come from HMG-box paralogs (example: VSP-0980, 30 reads at 47% identity to MAT1-2-1, 26 real MAT1-1-1 reads at 54%). Real hits at T4 lie near 54% identity, so any filter above 55% removes them.
- So with only distant references, blastx finds MAT1-2-1 but not MAT1-1-1 in this genus, and cannot be used for a call on its own.

VSP-0931 (the "Neither" strain of the samtools table; k-mer call `low_depth`): blastx calls MAT1-2 at T1 and T2 (28 and 73 reads) at only 60 to 66% identity to MAT1-2-1, while real *F. oxysporum* strains hit at 84 to 98%. Control depth is normal. The Kraken report assigns 69% of its reads to "Fusarium sp. MBC 181" and none to F. oxysporum among the top species. Reading: VSP-0931 is probably not *F. oxysporum*, and carries a divergent MAT1-2-1 that nucleotide matching misses. Not verified; the species is not established and the locus was not assembled to the end.

## Update 3: recruit reads, assemble, run `detect` (curator idea)
Method (`recruit.slurm`, job 29551936): blastx of all R1 and R2 reads (not capped) against the T2 + T4 MAT and flank proteins (81 proteins after removing duplicate headers, which also drops a few proteins that share species and gene name); keep both mates of every hit read; SPAdes `--careful`; `detect --taxid 5506` (frozen `run-94c1a3b` database, which contains F. oxysporum records) on the scaffolds.
Pilot on 4 strains (`results/2026-10-06_fola_reads_blastx/recruit/`):
| Strain | pairs recruited | assembly | `detect` | Right? |
|---|---|---|---|---|
| VSP-0980 (MAT1-1) | 1,479 | 4 scaffolds, 11.4 kb | MAT1-1 high, 4.7-kb contig | yes |
| 50a (MAT1-2) | 1,786 | 5 scaffolds, 11.3 kb | MAT1-2 high, 4.4-kb contig (also a 1.2-kb contig) | idiomorph yes; locus is 4.4 kb here against 12.1 kb with flanks in the genome assembly |
| VSP-0947 (both) | 1,072 | 20 scaffolds, 13.8 kb | MAT1-2 high on a 580-bp contig; MAT1-1 contigs (1.5 to 3.4 kb) withheld | no: "both" missed, MAT1-2 call on a 580-bp contig |
| VSP-0931 (divergent) | 1,262 | 6 scaffolds, 9.6 kb | none called; all 7 loci withheld | no call |
- The recruited assemblies cover the genes and a little around them, not the flanks. `detect` withholds short contigs without flanks (`modelled_gene_bar`). It works when a contig reaches the flanks (VSP-0980, 50a) and fails when it does not.
- Four strains is a pilot, not an accuracy estimate.
- Next test (not run): a second round, mapping the reads to the first-round contigs to extend into the flanks, then a new assembly.

## Update 4: simulated mixed samples (k-mer path, Fola panel)
Method (`scripts/simulate_mixes.sh`, `results/2026-10-06_fola_reads_mixes/`): 8,000,000 R1 reads per mixture; a MAT1-1 strain and a MAT1-2 strain are combined as the first n reads of each (process substitution, no temporary FASTQ). Minor and major roles are both tested (0 to 100% MAT1-2 reads). Three pairs, each of two strains with shared k-mer depth of 15 or more, chosen by depth (the best-covered strains, not a random draw):
VSP-1150 + VSP-0798, VSP-0777 + VSP-2032, VSP-1057 + JCP360. The first run stalled on one node; the missing fractions were re-run (jobs 29551712, 29552687). Table: `mixes_all.tsv`.

| MAT1-2 reads | VSP-1150+VSP-0798 | VSP-0777+VSP-2032 | VSP-1057+JCP360 |
|---|---|---|---|
| 0% | MAT1-1 (trace MAT1-2) | MAT1-1 (trace MAT1-2) | MAT1-1 (trace MAT1-2) |
| 1% / 2% | MAT1-1 (trace MAT1-2) | MAT1-1 (trace MAT1-2) | MAT1-1 (trace MAT1-2) |
| 5% | MAT1-1 (trace MAT1-2) | both (ratio < 0.25) | both (ratio < 0.25) |
| 10% | both (ratio < 0.25) | both (ratio < 0.25) | both (ratio < 0.25) |
| 20% / 50% | both | both | both |
| 80% | both | both (ratio < 0.25) | both |
| 90% to 99% | MAT1-2 (trace MAT1-1) | not run | MAT1-2 (trace MAT1-1) |
| 100% | MAT1-2 | not run | MAT1-2 |

- Pair 2 stops at 80%: VSP-2032 has about 7.0M R1 reads (fastp), so mixtures above 87% would not hold 8M reads and the true fraction would be wrong.
- Detection limit for `both` in these 3 pairs: a minor MAT1-2 share of 10% is called `both` in 3 of 3 pairs, 5% in 2 of 3, 2% and below in 0 of 3. A minor MAT1-1 share of 20% is called `both` in 3 of 3 pairs; 10% and below in 0 of 2 pairs (pair 2 was not run above 80%); those are called MAT1-2 with a `trace_MAT1-1` flag.
- The `trace_MAT1-2` flag is background: it appears at 0% MAT1-2 in all three pure MAT1-1 strains (MAT1-2 depth 1.0 to 1.6 against 15 to 17). Some part of the MAT1-2 reference occurs in these MAT1-1 genomes. The flag therefore does not detect a minor MAT1-2 idiomorph. `trace_MAT1-1` is not background: it is absent at 100% MAT1-2 and present at 90 to 99% (breadth 0.21 at 99% in pair 1).
- So a minor MAT1-1 share of 1 to 10% is visible only as `trace_MAT1-1`, and a minor MAT1-2 share below 5% is not visible. The `both` rule (relative depth 0.10) is set above this background; no change was made to it.
- Limits: three pairs; R1 only; fractions are shares of reads, not of nuclei; two haploid strains mixed, not a heterokaryon; no sequencing-error or index-hopping model. The A. fumigatus putative hybrids (Lofgren et al.) are not run.

## Update 5: second and later rounds of recruitment for Update 3 (curator request)
Method (`recruit_ext.slurm`, job 29600198): start from the round-1 scaffolds of Update 3. Each round maps all read pairs (R1 and R2, uncapped) to the previous scaffolds with minimap2 `-ax sr`, keeps pairs with a mate mapped at MAPQ >= 10, reassembles them with SPAdes `--careful`, and runs `detect --taxid 5506`.
Guards: stop after round 5, if recruited pairs grow less than 5%, or if more than 1,000,000 pairs. No run hit a guard; all four reached the round limit before converging (growth was still 12 to 18% per round except VSP-0947, below). Per-round tables: `results/2026-10-06_fola_reads_blastx/recruit_ext/<strain>/rounds.tsv`.

| Strain | Round 1 (Update 3) | Round 5 |
|---|---|---|
| 50a (MAT1-2) | 1,786 pairs; longest contig 4.5 kb; locus 4.1 kb | 6,035 pairs; longest 15.8 kb; MAT1-2 high, 15.5-kb locus with SLA2, APN2, COX13 |
| VSP-0980 (MAT1-1) | 1,479 pairs; 4.7 kb | 2,613 pairs; 15.3 kb; MAT1-1 high, 15.1-kb locus with the 3 flank genes |
| VSP-0931 (divergent) | 1,262 pairs; no call | 3,562 pairs; 15.6 kb; MAT1-2 high from round 2 (medium in round 3); 15.4-kb locus with the 3 flank genes |
| VSP-0947 (both) | 1,072 pairs; MAT1-2 only, 580-bp contig | 1,742 pairs; 5.8 kb; MAT1-1 high from round 2; MAT1-2 only in round 3 (1.2-kb contig, no flanks); no MAT1-2 in rounds 4 and 5 |

Checks against the full-genome assemblies (BLASTN, evalue 1e-50):
- 50a: the round-5 contig of 15,778 bp is 100.000% identical over its full length to NODE_86 of the independent whole-genome assembly (positions 115,322 to 131,099). The extension is correct.
- VSP-0980: the round-5 contig of 15,270 bp matches NODE_12 of the AVITI assembly (298,998 to 313,026) at 99.1 to 99.2%, which is 130 or so differences, not 0. I do not know the cause. The reads of VSP-0980 in `05_Pre-Processing` may not be the library of that assembly. This fits the lower MAT1-1 k-mer breadth of VSP-0980 (0.73, median 0.975 for MAT1-1 strains) when its own assembly supplied the panel. Not verified; it needs N. L.
- VSP-0931: no whole-genome assembly exists, so the 15.4-kb locus is not checked.

False second calls: round 5 also reports a second "high" MAT1-2 locus in 50a (NODE_3, 5.3 kb, 100% identical to another place in the genome, NODE_41) and in VSP-0980 (NODE_2, 2.3 kb, 98.7% to a non-MAT region of NODE_2 of the genome assembly). The VSP-0947 round-3 MAT1-2 call (1.2-kb contig) is the same type. In every real locus the three flank genes (SLA2, APN2, COX13) are present on the contig (3 of 3 loci in the table above); in these false calls they are absent or only one is present (3 of 3). That rule was read off six loci after the fact; it has not been tested on other strains.

What it shows:
- With extension, the recruit-and-assemble path reaches a full flanked locus in 3 of 4 strains, including VSP-0931 (no call before), and with no species-specific reference needed beyond the blastx proteins.
- It does not recover the second idiomorph of VSP-0947 (heterokaryon or mixed library): the MAT1-1 locus stays at 4.9 to 5.8 kb, and the MAT1-2 contig stays short. Mixed or two-allele samples are a failure case here.
- Without a flank-gene requirement, short paralog contigs are called as extra loci.
- Cost per strain: four rounds of minimap2 on the full read set plus four small SPAdes runs (8 CPUs). The run time per strain was not recorded.
- Limits: four strains, chosen as known cases (not a random sample); the extension rounds have no truth for VSP-0931 and VSP-0947; rounds stopped at the limit, not at convergence.

## Update 6: Aspergillus fumigatus, 304 assemblies and 331 read sets (assembly profile versus reads)
Question (curator): how well are MAT1-1 and MAT1-2 retrieved from assemblies (`detect`) and from raw reads (`reads-type`) in a large population set, and how do the two compare?
Branch `afum-reads-type` (stacked on `reads-type`, PR #45).

Data (`results/2026-10-07_afum_reads/strain_map.tsv`):
- Assemblies: `/bigdata/stajichlab/shared/projects/Afumigatus_pangenome/scaffolded/genomes/*.sorted.fasta`, 304 files.
- Reads: `/bigdata/stajichlab/shared/projects/Afum_popgenome/asm/input/`, 331 read sets from `asm/samples.dat` (column 1 = read-file prefix, column 2 = strain = genome name). 311 are paired and 20 are single-end (`_R1` only; typed as one file).
- 303 genomes have reads (Af10 has none); 28 read sets have no genome. 20 DMC and DMC2 pairs (for example DMC_AF100-1_3 and DMC2_AF100-1_3) give identical read results, so they are duplicates under two names.

Panel (`panels/aspergillus_fumigatus/`): the curated database records `746128_a1163_MAT_MAT1-1` (A1163, DS499596.1, 13,260 bp) and `746128_af293_MAT_MAT1-2` (Af293, NC_007196.1, 13,096 bp), whole locus segments with the flank genes. The two segments align at 99.3% over the flanks. The idiomorph-specific regions are MAT1-1 positions 4,576 to 6,664 (2,089 bp) and MAT1-2 positions 6,225 to 8,525 (2,301 bp) (`specific_regions.fasta`).

Assembly truth (`scripts/assembly_idiomorph_blast.py`, `blast/assembly_truth.tsv`): BLASTN of the two specific regions against each assembly (identity >= 90%, e-value 1e-10); present if >= 0.80 of the region is covered, absent if <= 0.20, otherwise partial. Result for 304 genomes: MAT1-1 only 160, MAT1-2 only 128, both 8, partial 7 (excluded from agreement counts), neither 1 (Moose_3-1.1).

**A panel flaw found and fixed.** With the first run (all unique k-mers), reads agreed with the assembly truth in 288 of 296 strains (97.3%), and 6 MAT1-2 assemblies were called `both`. In MAT1-2 strains the MAT1-1 "unique" k-mers had median breadth 0.351 and about half the MAT1-2 depth, as high as the signal in those 6 strains. Cause: the A1163 and Af293 flanks differ by SNPs, and every k-mer that spans a SNP counted as idiomorph-unique, and matches any strain with that flank allele. Fix (`Panel.from_sequences(min_run=100)`, option `--min-unique-run`, default 100): keep a unique k-mer only if it lies in a run of at least 100 consecutive unique positions (a SNP gives a run of at most 31). It removes 1,175 of 3,739 MAT1-1 and 1,168 of 3,575 MAT1-2 unique k-mers in this panel (31%), and 111 of 4,828 and 92 of 4,711 in the Fola panel (2%). Tests: 3 new, 47 pass with the smoke test. The rule was chosen after seeing the first result, so the 99.0% below is not an independent test of it. Fola is unchanged: 147/148 (`v4_concordance.tsv`).

Results with the filter (`comparison.tsv`; excluding the 7 partial and the 28 reads-only strains):
| Method | Agree with assembly truth |
|---|---|
| reads-type, 8M reads | 293/296 (99.0%); without the 20 DMC2 duplicates 280/283 |
| `detect` on assemblies (frozen run-2d1860a) | 288/297 (97.0%); without the DMC2 duplicates 275/284 |
- reads-type: MAT1-1 160/160, MAT1-2 126/127, both 7/8. Differences: AF100-12_2 (below), F18149-Manchester (below), Moose_3-1.1 (`low_depth`).
- `detect`: all 288 assemblies with one idiomorph got the right idiomorph (160 MAT1-1, 128 MAT1-2), 0 missed. All 8 assemblies with both regions got a single call (5 MAT1-1, 3 MAT1-2); the second idiomorph is not reported. All 304 genomes had at least one locus; 185 loci high and 120 medium confidence; no run failed.
- All reads (331): MAT1-1 174, MAT1-2 135, both 18, low_depth 4. Reads-only strains (28): MAT1-1 11, MAT1-2 7, both 7, low_depth 3.

The 8 assemblies with both regions (read depths, MAT1-1 / MAT1-2, k-mer units): 08-36-03-25 4.3/14.0, AF100-12_5 19.3/3.8, Afu_343_P_11 15.4/3.3, B7586_CDC-30 26.9/4.2, DMC_AF100-1_18 22.2/3.8, IFM_59359 24.0/16.4, IFM_61407 15.8/14.0: all `both` from reads. F18149-Manchester: reads call MAT1-1 (16.5) with MAT1-2 at 1.4 (breadth 0.71, flagged `trace_MAT1-2`), so the reads support the second region only weakly although the assembly holds both. Two of the three Lofgren et al. putative hybrids (IFM_59359, IFM_61407) are `both` from reads and from the assembly regions; `detect` reports MAT1-2 only (medium) for both. The third (AF100-1_3) has no assembly in this folder; its reads (DMC_AF100-1_3, DMC2_AF100-1_3) are `both`, at 13.3 / 13.4.

Reads show a second idiomorph that the assembly lacks: AF100-12_2 (assembly MAT1-2 only; reads MAT1-1 depth 3.5 against MAT1-2 18.8, MAT1-1 breadth 0.546). It is marginal: in the 126 MAT1-2 strains called MAT1-2 the MAT1-1 breadth has median 0.144 and maximum 0.443, and all 126 carry a `trace_MAT1-1` flag. In the 160 MAT1-1 strains the MAT1-2 breadth has median 0.000, maximum 0.514, and 46 carry a flag.

`detect` findings that are not read-typing results:
- 168 of the 169 MAT1-1 loci are labelled `homothallic_candidate`, and none of the 135 MAT1-2 loci are (`mat_locus`). The A1163 database record carries a degenerate MAT1-2-1 remnant as a core gene, and the class rule (`idiomorph.py`, both idiomorphs' core genes on one locus) fires on it. For *A. fumigatus* this label is therefore not evidence of homothallism. I did not change the record or the rule.
- Moose_3-1.1: `detect` reports MAT1-2 (high, scaffold_3:476,818-485,874) and MAT1-1 (medium, `partial_locus`, scaffold_3:1,319,318-1,326,086), 840 kb apart. Neither specific region matches the assembly at 90% identity, and its reads have shared depth 0.02, so the reads are not from an *A. fumigatus*-like genome. The species of this sample was not checked.

Limits: the assembly truth uses regions from the same two reference strains as the panel; it is independent of the reads but not of the reference choice. 20 duplicate read-set pairs inflate n. Only the first 8M reads were used; single-end samples (20) have half the reads of a pair at the same cap. A `both` call does not separate a heterokaryon, a diploid, a mixed culture or contamination. Reads-only strains have no truth. Cut-offs (breadth 0.50, relative depth 0.10, `min_run` 100) are untested outside Fola and this set.

Curator decisions needed: (1) keep `--min-unique-run` at 100. (2) The `homothallic_candidate` label on every *A. fumigatus* MAT1-1 locus: change the record (drop the remnant from the core genes), the rule, or leave and note. (3) Species of Moose_3-1.1. (4) Whether the second idiomorph in AF100-12_2, AF100-12_5, F18149-Manchester is a mixed sample (the AF100-12 series suggests one lab series) and who holds the cultures. (5) Whether `detect` should get the report-only "unlinked second idiomorph" check noted in the spec.

## Curator decisions
Open: (1) keep the 0.50 breadth fraction and 0.10 relative-depth cut-offs (now checked against simulated mixes, Update 4). (2) Alignment path: build, or accept k-mers only for same-species panels. (3) Whether to follow up VSP-0947 as a real two-idiomorph strain (heterokaryon or diploid) or a contaminated library: depths are 7.7 and 5.7 against shared 15.5, with assembled contigs both near 8x. (4) blastx for populations without a species locus: use T2-level references (same genus) with an identity filter of 40 to 80% (policy: a high filter removes paralog noise but loses divergent real alleles like VSP-0931); do not use distant references alone (T4). (5) Whether to assemble VSP-0931 to confirm it is a divergent MAT1-2-1 in a non-oxysporum Fusarium. (6) The second round of recruitment was built and tested (Update 5): RULED 2026-10-07 (J. Stajich): a contig needs at least 1 flank gene to be called a locus, probably 2; the exact number is not fixed and the rule is not implemented. More rounds and more strains are left for later; the pilot is accepted as a good start and the work moves on. (7) VSP-0980: ask N. L. whether the reads in `05_Pre-Processing` belong to the AVITI assembly used for the panel.

## Files
- `results/2026-10-06_fola_reads_type/` (samples.tsv, run_v2.slurm, v2_concordance.tsv); full per-strain TSVs in the main checkout `out_v2/`.
- `results/2026-10-06_fola_detect_all/`, `results/2026-10-06_fola_50a/` (reports, GFF3, scripts).
- `scripts/compare_reads_type.py`; code `src/MATPredict/reads/`.
- Assembly of 50a (not tracked): `/bigdata/stajichlab/jstajich/projects/fola_reads/50a_spades/`.
