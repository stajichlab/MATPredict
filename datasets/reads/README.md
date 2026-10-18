# Raw-read datasets for MATPredict reads typing

One folder per genus. Each holds a `manifest.tsv` (every read set found), notes on it, and, for the SRA genera, `selected.tsv` and `download.sh`.

| Folder | Source | Status |
|---|---|---|
| `coccidioides/` | Reads already on hpcc (`2025_All_Cocci/Assembly/input`), joined to the assembly freeze | manifest built; no reads copied |
| `cryptococcus/` | ENA/SRA, Illumina paired WGS | manifest under review |
| `fusarium/` | ENA/SRA, Illumina paired WGS (*F. graminearum* complex and other species) | manifest under review |

## Coccidioides

`build_manifest.py` joins three hpcc files (read-only copies; the commands are in the script header):

- `.../2025_All_Cocci/Assembly/samples.csv` - one row per read set: strain, centre, locus tag, and the species label of that file
- `.../NovInvenio_Investigations/studies/fungi/coccidioides_pangenome/species.csv` - corrected species and the freeze assembly path for 529 strains
- `ls -L` of `.../Assembly/input` - file names and sizes

Species corrections (`coccidioides_pangenome/CORRECTIONS.md`, 2026-09-21), shown with `relabelled = yes`:
`485B-1_L_OLD_CPA0023` immitis to posadasii, `B3476` posadasii to immitis. `samples.csv` still carries the old labels and the freeze file names (`Coccidioides_<species>_<strain>.*`) also still carry the wrong species word; `species.csv` has the corrected labels.

Counts: 558 read sets (383 *C. posadasii*, 175 *C. immitis*), about 1.23 TB of gzipped FASTQ. 529 have a row in the assembly freeze (`species.csv`); 29 do not.
The freeze's `DNA/*.scaffolds.fa` links all dangle (the target `Assembly/unscaffolded_annotation/` no longer exists), so `manifest.tsv` takes `assembly` from `Assembly/asm/AAFTF/<strain>.sorted.fasta` instead: 503 of 558 strains have one, and the dangling freeze path is kept in `freeze_assembly_link`. This is the unscaffolded assembly that fed the annotation; I did not check that it is byte-identical to what the freeze once pointed at.
`TX13` is labelled *Fusarium solani* in `samples.csv`; it is dropped (decided 2026-10-17).

## SRA datasets: cap and data handling

- Cap: 50 runs per species (decided 2026-10-17). `fusarium/selected.tsv` is the rank <= 50 set (960 runs); the earlier 150-cap set is `selected_extended.tsv`. `cryptococcus/selected.tsv` (312 runs) was already capped at 50.
- Reads are not kept on /bigdata. `fetch_type.slurm` downloads one run into `$SCRATCH`, types it against a species panel, keeps only the result TSV and deletes the reads. It is syntax-checked only, and panels for Cryptococcus and Fusarium do not exist yet (Coccidioides panels are in `panels/coccidioides/`).

## Coccidioides pilot result (24 strains, 12 per species)

Per strain, `matpredict detect` on the assembly versus `reads-type` (first 8M reads) on each panel; table in `coccidioides/pilot/comparison.tsv`.

- `detect` called every strain: 12 MAT1-1 and 12 MAT1-2, 22 of 24 with high confidence.
- immitis panel (`panels/coccidioides/immitis_*`): 23 of 24 reads calls agree with `detect`, in both species (12 of 12 *C. posadasii*, 11 of 12 *C. immitis*). The exception is CA10, `low_depth` with 0.46x MAT1-1 depth.
- posadasii panel: 0 of 24 called, all `low_depth`. Its MAT1-1 (EF512013) and MAT1-2 (EF472258) records are separate deposits with no shared flank, so the shared-depth control is 0. The idiomorph depths are right (they match `detect`); the panel needs a shared region added before it is usable.
- The two relabelled strains behave as their corrected species: both are called correctly by the immitis panel, B3476 (relabelled *immitis*) is typed MAT1-1 and 485B-1_L_OLD_CPA0023 (relabelled *posadasii*) MAT1-2, in agreement with `detect`.
- Noise: every strain also carries a `trace_` flag for the other idiomorph, often at a substantial depth (e.g. MAT1-1 55x vs MAT1-2 42x in B11343, a MAT1-2 strain), because the immitis panel includes flanking sequence that is not identical between strains. The calls are right, but the panel should be trimmed to the idiomorph-specific region.

Note on subsampling: old SRA FASTQs (2012-13 HiSeq) are position-sorted, so the first N reads can all come from one chromosome. `fetch_type.slurm` therefore keeps every Nth read pair (about 8M pairs, N from the manifest's `read_count`) rather than using `--max-reads`. The Coccidioides pilot used the first 8M reads of unsorted centre files and was fine.

### Trimmed panel (`panels/coccidioides/cocci2_*`), same 24 strains

`coccidioides/pilot/comparison_cocci2.tsv`: 23 of 24 reads calls agree with `detect` (CA10 is the exception, `low_depth`, with 0x depth on every panel, including the shared flank, so its read set does not carry the locus at all; possible read/assembly mix-up, not checked). The wrong idiomorph now reads 0.00-0.15x in every strain (it was up to 100% of the true depth with the untrimmed panel).

## Coccidioides full run (558 strains, `cocci2` panel)

Table: `coccidioides/results/cocci2_reads_vs_detect.tsv` (reads call, triage, depths, and the `detect` call and confidence for the 503 strains with an assembly). Reads: first 8M reads per strain; `detect` v0.6.x main (21d7bfb) on the AAFTF assemblies. All 1,116 jobs completed.

| reads call | *C. immitis* | *C. posadasii* |
|---|---|---|
| MAT1-1 | 78 | 189 |
| MAT1-2 | 95 | 184 |
| both | 1 | 3 |
| low_depth (omit) | 1 | 7 |

- 546 strains typed to one idiomorph, 4 `both`, 8 `low_depth`. The `low_depth` strains have little or no shared-flank depth (0.0 to 7.3x) and are omitted as low coverage (CA10 is one).
- Of the 495 typed strains that also have a `detect` result, 487 agree. The 8 that differ: 6 where `detect` found no locus in the assembly (the reads call MAT1-1 or MAT1-2 at shared depth 1.1 to 4.0x), 2 where `detect` called both idiomorphs in the assembly while the reads show one (NM_0071, TX7). Not traced.
- The 4 `both` reads calls are not checked against a second source. Whether they are mixed samples, heterokaryons or panel cross-talk is open.
- 55 strains have no assembly, so they have a reads call only.
