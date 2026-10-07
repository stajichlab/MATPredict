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
- Mixed samples are not calibrated yet. The planned simulated mixes (50:50, 80:20, 95:5) are not run. A 2% contaminant (like 50a's trace) would fall under the 0.10 depth cut-off and be flagged, not called.
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

## Curator decisions
Open: (1) keep the 0.50 breadth fraction and 0.10 relative-depth cut-offs, or recalibrate on simulated mixes first. (2) Alignment path: build, or accept k-mers only for same-species panels. (3) Whether to follow up VSP-0947 as a real two-idiomorph strain (heterokaryon or diploid) or a contaminated library: depths are 7.7 and 5.7 against shared 15.5, with assembled contigs both near 8x.

## Files
- `results/2026-10-06_fola_reads_type/` (samples.tsv, run_v2.slurm, v2_concordance.tsv); full per-strain TSVs in the main checkout `out_v2/`.
- `results/2026-10-06_fola_detect_all/`, `results/2026-10-06_fola_50a/` (reports, GFF3, scripts).
- `scripts/compare_reads_type.py`; code `src/MATPredict/reads/`.
- Assembly of 50a (not tracked): `/bigdata/stajichlab/jstajich/projects/fola_reads/50a_spades/`.
