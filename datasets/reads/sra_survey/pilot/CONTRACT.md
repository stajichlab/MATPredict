# Interface of the raw-read typing pipeline (runner-independent)

The unit of work is `genotype_mat_one_accession.sh`: one run in, two small files out. Slurm arrays (`stream_type_sra.slurm`, `submit_species.sh`) and a Nextflow process call the same script, so changing runner changes nothing in the science. Nothing in the worker reads Slurm variables.

## Input: the run list (`runs.tsv`, tab-separated, header row)
| Column | Required | Meaning |
|---|---|---|
| `run_accession` | yes | SRR, ERR or DRR accession; one row per strain (the survey picks the run) |
| `read_count` | yes | number of spots (ENA `read_count`); the windows are laid out along it. 0 or blank: read the whole run |
| `panel_dir`, `panel_prefix` | yes (here, or as arguments) | `PANEL_DIR/PANEL_PREFIX_<IDIOMORPH>.fasta` per idiomorph; a row may name its own panel so one list can mix species |
| `group` | no | free label (species, dataset); the summary groups by it |
| `library_layout` | no | PAIRED or SINGLE; recorded only |
| `strain`, `study_accession` | no | carried into the summary |
| `truth_detect`, `truth_assembly` | no | an independent idiomorph call (the `detect` call of an assembly with the same BioSample) and the assembly it came from; used for agreement |

## Worker: `genotype_mat_one_accession.sh RUN SPOTS PANEL_DIR PANEL_PREFIX OUT_DIR [GROUP] [LAYOUT]`
Environment: `TARGET_SPOTS` (2,000,000 spots), `WINDOWS` (8), `SRA_URL`, `IDIOMORPHS` ("MAT1-1 MAT1-2"), `SCRATCH` (temporary space), `KEEP_SRA`, and the tool locations (`SRATOOLKIT_BIN`, `PY`, `PYTHONPATH`, `MATPREDICT_DB_ROOT`). Needs `curl`, the SRA toolkit and `python -m MATPredict reads-type`. Network: HTTPS to `sra-pub-run-odp.s3.amazonaws.com`.

- **Output** in `OUT_DIR`: `<run>.<prefix>.tsv` (the `reads-type` row: `call`, `MAT1-1_breadth`, `MAT1-2_breadth`, `MAT1-1_depth`, `MAT1-2_depth`, `shared_depth`, `reads_used`, `flags`) and `<run>.run.json` (`run`, `group`, `layout`, `spots_total`, `target_spots`, `windows`, `sra_bytes`, `secs_download`, `secs_convert`, `secs_type`, `state`, `source`, `host`).
- **State** in the json: `done`, `unavailable_at_ncbi` (the run is in ENA but neither NCBI's AWS copy nor NCBI itself has it; about 1 run in 100 so far), `failed_download`, `failed_convert`, `failed_type` (`started` if the task died before writing a state). Calls: `MAT1-1`, `MAT1-2`, `both`, `none`, `low_depth` (idiomorph names follow the panel).
- **Exit**: 0 done or already done, 1 failed, 2 bad arguments. **Idempotent**: the `.tsv` is the done-marker, so any runner may re-run the whole list and only the missing runs work. **No shared state** between runs; each task has its own temporary directory, removed on exit. The reads are never kept.
- **Resources** (measured on 20 runs): about 1 GB download, 110 s, 5 GB peak memory (ask for 8), up to 4 CPUs in the typing step, a few GB of temporary space.

## Summary: `summarize_pilot.py --runs RUNS.tsv --results OUT_DIR --out SUMMARY_DIR`
Reads only the small files above; rerun any time. Writes `SUMMARY.md`, `summary_by_group.tsv`, `results_all.tsv`, `discordant.tsv`, `failures.tsv`.

## Runners
- **Slurm** (now): `submit_species.sh SLUG RUNS.tsv` queues an array, a retry array and a summary job; results and logs under `ROOT/SLUG/`.
- **Nextflow** (kept possible, not used): `main.nf.example` shows the same list as a channel and the worker as one process (`maxForks` for the concurrency limit, `errorStrategy retry`, `publishDir` for the outputs). It has not been run.
