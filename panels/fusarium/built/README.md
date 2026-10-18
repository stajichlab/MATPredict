# Missing-idiomorph candidate panels for Fusarium (built from public reads, 2026-10-17/18)

For species whose public assemblies show one idiomorph only, the missing idiomorph was recovered from ENA reads of ONE strain
(recruit, assemble, extend), then cut together with its shared flank into a two-record panel in the repo's format
(`<species>_MAT1-1.fasta`, `<species>_MAT1-2.fasta`, one record each). Use as in `../README.md`:

```
matpredict reads-type --idiomorph MAT1-1=panels/fusarium/built/andiyazi_MAT1-1.fasta \
                      --idiomorph MAT1-2=panels/fusarium/built/andiyazi_MAT1-2.fasta --reads R1.fq.gz R2.fq.gz --out out.tsv
```

Files: `validation.tsv` (every read test, 68 rows), `screen_results.tsv` (which idiomorph each screened run carries, from reads, independent of any panel),
`species_missing_idiomorph.tsv` (step 1 list), `annotation_tblastn.txt` (MAT proteins found in each record), `scripts/` (every script run on hpcc, plus a copy of the repo's
`reads` typing code as used), `round*.list` (job inputs).

## Result in one table

| species | MAT1-1 record (source) | MAT1-2 record (source) | status |
|---|---|---|---|
| F. andiyazi | `GCA_044706615.1:JADOXN010000013.1:364752-370443` (5692 bp, assembly YN28) | `ERR4080120:NODE_1_length_13673_cov_20.229182:5025-9929(-)` (4905 bp, FFSC-wgs08) | validated on 9 of 10 other-run tests; 1 run no call (see below) |
| F. keratoplasticum | `GCF_025433545.1:NC_070531.1:3466907-3477594` (10688 bp, Fu6.1) | `SRR37004787:NODE_1_length_13853_cov_36.002186:4052-10436` (6385 bp, B27623) | 8/8 correct, but 5 of the 6 MAT1-2 runs are one Israeli clinical project and only 2 MAT1-1 runs exist |
| F. avenaceum | `SRR19681895:NODE_1_length_8761_cov_23.440120:1402-6603` (5202 bp, NRRL 36457) | `SRR29914801:NODE_1_length_6214_cov_26.899821:538-4862` (4325 bp, INRA874) | 8/8 correct (only 2 MAT1-2 runs: source + SRR29914798) |
| F. brachygibbosum | `SRR37004795:NODE_1_length_13724_cov_53.432448:4035-10248` (6214 bp, B27615) | `SRR37004789:NODE_1_length_13170_cov_55.243042:3842-9325(-)` (5484 bp, B27576) | NOT validated: 1 of 2 independent runs typed (B27584 no call, B27603 call at breadth 0.50) |

"Validated" here means: other runs of the species (not the source run) typed to the idiomorph that the reads themselves carry. Only runs outside the source
strain count; source/same-strain runs are listed in `validation.tsv` (role column) but are circular.
Both records of a panel come from different strains of the species, except brachygibbosum and avenaceum where both come from reads (different strains each);
the MAT1-1 records of andiyazi and keratoplasticum come from the NCBI assembly of the species.

Coordinates are 1-based inclusive. A `(-)` means the record is the reverse complement of the named contig (so both records of a panel read in the same direction).
Contig names are SPAdes names from the recruited-read assemblies (`NODE_n_length_L_cov_C`, C = SPAdes k-mer coverage); contigs are not deposited.

## Species list (step 1)

`species_missing_idiomorph.tsv`: from `species_panel_map.tsv`, 46 species (186 selected runs) have one assembled idiomorph and 37 species (100 runs) have no assembly.
Those with at least 3 selected runs and one idiomorph: andiyazi 36, venenatum 27, avenaceum 21, keratoplasticum 17, euwallaceae 8, phialophorum 6 (has the fola panel),
kuroshium 6, nygamai 5, brachygibbosum 4, virguliforme 4, commune 3. No-assembly species (Fusarium sp. 31, petroliphilum 13, duoseptatum 5, ambrosium 4, cf. hostae 3, nematophilum 3,
grosmichelii 3, tucumaniae 3, ...) were not touched. Not attempted for lack of budget/large downloads: virguliforme (4 runs of 64-77 M pairs), commune, phialophorum, and every species with 1-2 runs.
Check of the earlier claim: avenaceum is NOT MAT1-1-only. Assembly GCA_025948275.1 (WV21P1A) carries a full MAT1-2 locus (CP109669.1:1220506-1224798, 96.7% identical to my read-built
MAT1-2 record; MAT1-2-1 is only 55-58% identical to the fujikuroi query, which is why it was missed). It is an independent confirmation of the avenaceum MAT1-2 record.

## What was screened and what was built

Every screen is `diamond blastx` of an every-Nth-pair subsample (about 2 M pairs) against the MAT proteins of the repo db (`db/Ascomycota/Hypocreales/*/proteins.faa`: F. fujikuroi MAT1-1-1/-2/-3 and MAT1-2-1,
F. graminearum MAT1-1-x/MAT1-2-1, and two more), counting reads per gene (`screen_results.tsv`). MAT1-2-1 hits at 30-41% identity are HMG-box noise.
Findings per species (reads decide, not assemblies):

| species | runs screened | MAT1-1 | MAT1-2 | outcome |
|---|--:|--:|--:|---|
| andiyazi | 10 | 5 | 5 | built |
| keratoplasticum | 4 screened + 4 more in validation | 2 (LHS11, Fu6) | 6 (5 Israel PRJNA1391476, 1 India) | built |
| avenaceum | 5 screened + 3 more in validation | 6 | 2 (INRA874, INRA497) | built |
| brachygibbosum | 2 + 2 | 2 (B27615, B27603) | 2 (B27576, B27584) | built, not validated |
| venenatum | 3 (all runs are strain A3/5) | 3 | 0 | cannot build: single strain, MAT1-1 only |
| nygamai | 4 of 5 | 4 | 0 | cannot build from the selected runs |
| euwallaceae | 7 of 8 | 7 | 0 | cannot build (MAT1-2-1 at ~45-60% to the query, detection is weak at this divergence) |
| kuroshium | 3 of 6 | 3 | 0 | cannot build from screened runs |

## Method and commands (hpcc, Slurm `short`, account `stajichlab`; all scripts in `scripts/`)

Tools (modules unless noted): diamond 2.1.12, NCBI BLAST+ 2.17.0+ (blastn, tblastn, makeblastdb), seqtk 1.4, aria2 1.36.0 (ENA https downloads; plain curl/ftp from the nodes ran at about 0.5 MB/s and
stalled, aria2 gave 5.8 MB/s), apptainer + **AAFTF image `AAFTF-0.7.0-beta.4-18-g9ce56a1.sif`** (`/bigdata/stajichlab/shared/singularity_cache/`, the `NXF_APPTAINER_CACHEDIR` path; `$NXF_APPTAINER_CACHE` is not set)
containing AAFTF 0.7.0-beta.4, SPAdes 4.3.0, bwa 0.7.19, samtools 1.24. Python 3.9 (miniconda module) for the typing code. Not used: megahit, minimap2, miniprot, exonerate (installed, not needed).

1. `prep.sh`: 6-frame translation (`orfdb.py`, ORF segments >= 40 aa) of every Fusarium panel record in `panels/fusarium/` plus the fola panel, joined with the repo MAT proteins, `diamond makedb` (recruit db); NCBI assemblies via `datasets` API.
2. `screen.sh RUN`: download with md5 check (ENA `fastq_ftp` md5 in `runs.tsv`), every-Nth read pair (N = read_count / 2,000,000, never the first N), `diamond blastx -e 1e-5 -k 1` against the MAT proteins.
3. `recassemble.sh RUN 3` (inside Slurm, 16 cpu, 48 GB): `diamond blastx -e 1e-5` of all reads (both mates) against the recruit db; `seqtk subseq` of both mates for the hit ids (so pairs stay together);
   `apptainer exec <sif> AAFTF assemble --method spades -c 16 -m 40` on the recruited pairs; contigs >= 300 bp with a blastx hit to the recruit db are the locus contigs; all reads are mapped to them
   (`apptainer exec <sif> bwa mem`, `samtools view -F 4 -q 20`, earlier runs `-q 1`), mapped pairs are added to the recruited set and the assembly is repeated, 3 extension rounds.
   Recruited pairs after round 0 and round 3: andiyazi 1280 to 2161; avenaceum MAT1-1 791 to 2071; brachygibbosum 1776 (MAT1-1) / 1642 (MAT1-2) to 3092 and 3060 by round 1;
   keratoplasticum 2383 to 41,863 (5 locus contigs, one real). **avenaceum MAT1-2 (SRR29914801) ran with `-q 1`: the recruited set exploded to 744,543 pairs at round 1 and 1,077,141 at round 3** (repeat contig); the final
   6214 bp locus contig is nevertheless confirmed against assembly WV21P1A (96.7% over 4323 bp) and by typing. The cap (stop above 300,000 pairs) and `-q 20` were added afterwards and used for the later runs.
4. Record cutting (`build_record.py`, same design as the repo's panels): `blastn` of the MAT1-1 locus against the MAT1-2 contig; the largest gap between the two long colinear blocks is the idiomorph core; each side gets up to 600 bp
   of flank where the two loci are >= 97% identical. Alignment columns next to the core are often diverged, so the core is extended to the first 100-column window >= 97% identical before the flank starts (without this,
   the right flank was dropped in all pilot cuts). Resulting flanks (left/right bp): andiyazi 600/450, keratoplasticum 600/600, avenaceum 325/600 (the MAT1-2 contig starts 49 bp before the flank block), brachygibbosum 600/600.
   Flank blocks align at 97.3-99.2% between the two source loci (andiyazi 98.4/98.5, keratoplasticum 98.9/98.5, avenaceum 97.3/98.6, brachygibbosum 99.2/98.1; blastn over the full blocks of 0.8-4.7 kb).
   The earlier version of the andiyazi and keratoplasticum records (right flank dropped) was typed once for andiyazi only, results not kept; all numbers below are for the final records.
5. Typing (`typejob.sh`): download, every-Nth pair to about 1.6 M pairs (3.2-3.9 M reads used), the repo's `reads-type` code unchanged (`scripts/typer/`, copy of `src/MATPredict/reads/{kmers,panel,typing}.py`, k=31, min-unique-run 100, run through
   `type_run.py` because the package is not installed on the cluster), plus the same-subsample diamond screen as independent truth. Reads are deleted at job end (job-local `$SCRATCH`); nothing but small result files is kept under
   `/bigdata/stajichlab/jstajich/projects/MATPredict_cocci/reads_asm/`.
6. Assembly check, andiyazi source (`covcheck.sh`, bwa in the AAFTF image): 1.97 M pairs mapped back to the 13,673 bp contig: mean depth 10.7x, minimum 1x, 0 uncovered bases, 1120 of 1174 mapped reads properly paired.
   This check was run for andiyazi only.

## Annotation evidence (`annotation_tblastn.txt`, `tblastn` of the repo MAT proteins on the final records)

- andiyazi MAT1-1: MAT1-1-1 (97.6% over 291 aa), MAT1-1-2 (97.5% over 401 aa), MAT1-1-3 (96% over 50 aa, then 84%). MAT1-2: MAT1-2-1 HMG-box and exons (98.5% over 66 aa, 96.6% over 89 aa, 91% over 58 aa of the F. fujikuroi MAT1-2-1).
  The MAT1-2 contig as a whole is 95.3% identical over 5035 bp to the F. verticillioides MAT1-2 record; flanks 98.5% to the YN28 MAT1-1 locus on both sides.
- keratoplasticum MAT1-1 (Fu6.1): MAT1-1-1/-2/-3 at 36-41% identity (solani complex is diverged from the fujikuroi queries). MAT1-2: MAT1-2-1 HMG hits 40-47% over 107-191 aa; records are 97.7/96.8% identical to the Fu6.1 flanks on both sides.
- avenaceum MAT1-1: MAT1-1-1 44-47% over 343-388 aa, MAT1-1-2 ~47%, MAT1-1-3 57%. MAT1-2: MAT1-2-1 55-60% over 40-100 aa pieces (HMG box and exons) plus the matching WV21P1A locus.
- brachygibbosum MAT1-1: MAT1-1-1 65% over 316 aa, MAT1-1-2 68% over 319 aa, MAT1-1-3 63% over 225 aa. MAT1-2: MAT1-2-1 55-75% over 67-201 aa.
- Flank genes (SLA2, APN1/APN2, COX13) were not annotated by name. Evidence that the contigs are the right locus: the MAT genes sit between two blocks that align at 97-99% to the other idiomorph's locus,
  and the contigs (13.2-13.9 kb; 8.8 and 6.2 kb for avenaceum) span the whole locus.
- MAT1-2-x accessory genes: not searched for (only MAT1-2-1 proteins exist in the repo db).

## Validation (`validation.tsv`)

Own panel, other runs of the species, calls from the repo code (1.6 M pairs each):

| species | independent MAT1-1 runs correct | independent MAT1-2 runs correct | no call / wrong |
|---|---|---|---|
| andiyazi | 3 of 4 (YN-SJ67, YN-SJ48, YN-DH30 MAT1-1, breadth 0.998-1.0, depth 7-9) | 4 of 4 (YN-SJ65, GM64, DH16, DH04: MAT1-2, breadth 0.71, depth 5-6) | FJ-CT46 (SRR22163438): `low_depth`, shared depth 0.57, MAT1-1 breadth 0.08 |
| keratoplasticum | 1 of 1 (LHS11, breadth 0.99, depth 6.6) | 5 of 5 (breadth 0.98-1.0, depth 7-13) | none |
| avenaceum | 5 of 5 (breadth 0.82-0.99, depth 4-17) | 1 of 1 (INRA497, breadth 0.98, depth 4.5) | none |
| brachygibbosum | 1 of 1 (B27603 `MAT1-1`, breadth only 0.50, shared depth 3.8) | 0 of 1 (B27584: `none`, MAT1-2 breadth 0.47 vs 0.50 needed at depth 5.1, `trace_MAT1-2`) | B27584 |

Source-run and same-strain tests (andiyazi ERR4080120, FA-YN28; keratoplasticum SRR37004787, Fu6 SRR18390558; avenaceum both sources; brachygibbosum both sources) all typed correctly, as expected.
Other-species panels (verticillioides for andiyazi; solani for keratoplasticum; fujikuroi, proliferatum, circinatum, xylarioides, sacchari for FJ-CT46) gave `none`/`low_depth` on every run tested: reuse across species did not work for these.

Cross-idiomorph k-mers: the `trace_` flags are NOT absent. The correct call carried a flag of the other idiomorph in andiyazi 3 of 9 runs (breadth 0.02-0.05, depth <= 0.05; the gene screen of those runs shows 1 read of the other MAT1-1 genes),
keratoplasticum in all 6 MAT1-2 runs (`trace_MAT1-1`, breadth 0.011-0.013, depth 0.06-0.2; the gene screen shows 0-1 reads of MAT1-1 genes in those runs), avenaceum in 2 of 8. These are below the call threshold and none changed a call.
Their origin (a repeated element in the solani MAT1-1 core, or barcode bleed) was not examined.

## Unverified / caveats

- Every panel has one source strain per idiomorph; allele variation inside a species is not represented. FJ-CT46 (andiyazi label; MAT1-1-1 96% to F. fujikuroi, MAT1-1-2 83%) is not typed by the andiyazi panel nor any fujikuroi-complex
  panel, so it is a divergent lineage or a mislabelled strain; B27584 and B27603 (brachygibbosum) show the same pattern at smaller scale. brachygibbosum is not validated; the label of the Israeli clinical set (PRJNA1391476, 13 keratoplasticum, 4 brachygibbosum, 1 nygamai) was not checked.
- Heterothallism was not tested beyond the screens: in every screened run only one idiomorph carried reads (no run called `both`), consistent with heterothallism, but the run counts are small and the keratoplasticum MAT1-2 runs are clonal-looking (one project).
- euwallaceae, kuroshium, nygamai, venenatum: MAT1-2 absence is only "not found in the screened runs"; for euwallaceae/kuroshium the MAT1-2-1 query identity is probably too low for blastx at 1e-5 on 2 M pairs.
- The avenaceum records have a 325 bp left flank, so their shared control is smaller than in the other panels (shared k-mer count not computed).
- MAT1-2-x genes, intron structure, and flank gene names were not annotated. Assembly stats are from the recruited-read assemblies only, not whole genomes; only andiyazi has a read-mapping check.
- Typing used 1.6 M pairs per run (3.2-3.9 M reads), fewer than the 8 M pairs used for other panels; breadths below 1 for MAT1-2 in andiyazi (0.71) are an allele effect or low depth, not examined.
- Unscreened runs of each built species were not typed (andiyazi 26 of 36, avenaceum 13 of 21, keratoplasticum 9 of 17 runs not typed).

## Compute

About 114 Slurm jobs on `short` (screens 61 incl. 24 failed downloads and 11 cancelled by me after a stalled-ftp problem, typing 42, assembly 6, other 5), roughly 250 CPU-hours (screens 161, assembly 60, typing 29),
downloads of an estimated 190 GB in total (not measured; 3-30 M pairs per run) that were deleted when each job ended. Round sizes were above the requested 20 jobs because job arrays were throttled (`%2-%6`) to avoid ENA connection limits.
