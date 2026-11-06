# Reads availability: what else could be MAT-typed from short reads

A refreshable join of **what read data exists** (public archives and the lab's own disks) against **what MATPredict has already typed**, per species. It answers: for which species can an idiomorph frequency be estimated from more isolates than we have used, and what blocks the rest?

Nothing is downloaded. Each refresh writes one dated snapshot under `snapshots/`.

## Sources

| What | From | Owner |
|---|---|---|
| Public runs (ENA `read_run`, Illumina WGS GENOMIC, all Fungi) | PopSeq `cache/fungi_runs_ena_<date>.tsv.gz` | PopSeq |
| Taxid to species (varieties and *f. sp.* rolled up to the species taxid, current NCBI name) | PopSeq `catalog/taxonomy_species_map.csv.gz` | PopSeq |
| Genome size per species (coverage estimate only) | PopSeq `catalog/genome_sizes.csv` | PopSeq |
| Local read sets under `Population_Genomics` | PopSeq `scripts/popseq_inventory.py` (spec 005), rerun at each refresh | PopSeq |
| Local subproject to species, and whether it is a population sample | `local_subprojects.tsv` (here, hand-curated) | MATPredict |
| Runs typed, calls | `sra_survey/pilot/wave1/results_all_final.tsv`, `wave2/results_all_final.tsv`, `pilot4/results_all.tsv` | MATPredict |
| Strains typed from local files | `coccidioides/results/cocci2_reads_vs_detect.tsv` | MATPredict |
| Wave gate per species | `sra_survey/pilot/wave{1,2}/gates_final.tsv` | MATPredict |
| Reference tier | `db/` (as `survey_ena.py`) | MATPredict |

PopSeq is the lab's shared catalog (`/bigdata/stajichlab/shared/projects/Fungi/PopSeq`, github `stajichlab/PopSeq`). It replaces MATPredict's own ENA query of 2026-10-26 (`../sra_survey/`): one harvest serves every project, and its taxid-keyed species fixes the synonym splits the name-keyed survey had.

## Rules

- **Isolate** = one BioSample. Per isolate the run nearest 3 Gb inside 0.6 to 12 Gb is chosen, else the largest (as `select_runs.py`).
- **Eligible** = the study is not flagged by `classify_study()` (the survey's `mutant`, `evolution`, `cross`, `not_wgs` rule, now shared from `survey_ena.py`), the run has at least 300 Mb, paired or single-end. *S. cerevisiae* is skipped (curator ruling 2026-10-26).
- **Typed** = a `done` run in a result table, mapped to its BioSample through the run table. Coccidioides strains typed from local files are matched to public isolates by strain name.
- **Local** = a sample of a counted subproject with at least one file. `local_likely_public` matches it to the archive by run accession, BioSample or strain name; `local_not_public` is therefore an **upper bound** on new strains (a strain-name match is weak, spec 005 section 4, and many local sets were deposited under other names).

## Action per species (`species_availability.tsv`, column `action`)

| action | meaning | what to do |
|---|---|---|
| `expand` | wave gate PASS or PROVISIONAL; untyped isolates exist (the 50-run cap) | type `untyped_runs.tsv.gz` rows for the species with its existing panel |
| `unvalidated` | typed, but the gate could not validate | curator ruling first (decision 1 in the reads-gate notes) |
| `held` | HOLD for typing rate or agreement | re-test or panel work |
| `typed_no_gate` | typed with a panel outside the waves (pilot4: *Pyricularia oryzae*) but never gated | gate it, then treat as `expand` |
| `blocked_panel` | HOLD because no panel could be built | reference or panel work, not more reads |
| `new_candidate` | not in a wave; records for the species or its genus | next wave (`wave_species.py`) |
| `needs_reference` | not in a wave; no record in the genus | reference curation priority, ranked by isolates |
| `homothallic` | passed as homothallic | none; no idiomorph frequency |
| `complete` | every eligible public isolate typed | none until the archive grows |

Also per species: eligible, typed and untyped counts; studies and `largest_study_frac` (near 1 = one study dominates, often clonal outbreak sampling, so not a population frequency); the current calls; and a binomial (Wilson) 95% half-width for the minor idiomorph now and if every untyped isolate were typed at the species' current typing rate (`ci95_halfwidth_now`, `ci95_halfwidth_all`), which ranks where more typing changes an estimate most.

## Files in a snapshot (`snapshots/<PopSeq run table date>/`)

- `SUMMARY.md`: counts per action and the top species of each.
- `species_availability.tsv`: one row per species.
- `local_samples.tsv`: every local sample PopSeq's inventory found, with its species and whether it is counted.
- `untyped_runs.tsv.gz`: one run per untyped eligible isolate of `expand`, `unvalidated`, `held` and `typed_no_gate` species; the run list columns of `../sra_survey/pilot/CONTRACT.md` once a panel is named. **Not committed** (about 1 MB, regenerated from the inputs).
- `snapshot.json`: input files with sha256, the PopSeq commit, the counts of runs kept and dropped.
- `CHANGES.md`: against the newest earlier snapshot: species new to the table, changed action, gained isolates, and **lost isolates** (a red flag: archives rarely drop records, so a loss usually means the query or a filter changed).
- `local_inventory_failed.txt`: subprojects PopSeq's inventory could not read.

## Refresh

On the HPCC, from the root of a MATPredict checkout:

```bash
sbatch datasets/reads/availability/refresh.sh     # or bash, about 2 minutes
```

It reruns PopSeq's local inventory (read-only on `Population_Genomics`), takes PopSeq's newest kingdom run table, and writes the snapshot with `CHANGES.md`. Review, then commit the snapshot folder. A subproject that PopSeq lists but `local_subprojects.tsv` does not map is printed and not counted: add it to the map.

**When**: after each PopSeq kingdom refresh (new `cache/fungi_runs_ena_<date>.tsv.gz`), and after each MATPredict typing wave is merged (the typed side changed). Both are events, so the refresh is run by hand at those points rather than on a timer. If a timer is wanted later, a monthly `scrontab` line is enough (`0 6 1 * * cd <checkout> && bash datasets/reads/availability/refresh.sh`), since the public archive grows by a few percent a month.

A new typing wave adds its `results_all_final.tsv` and `gates_final.tsv` to `DEFAULT_TYPED` and `DEFAULT_GATES` in `reads_availability.py` (or `--typed` / `--gates`).

## Known gaps (2026-11-06)

- PopSeq's kingdom run table sits in its untracked `cache/`; the snapshot records its sha256 but the file itself is not versioned by PopSeq. Asked of PopSeq: keep the dated table, or track it.
- PopSeq's inventory fails on dangling symlinks in `Coccidioides/C_immitis_WAref`, `Friedmanniomyces_endolithicus` and `Rhizopus_stolonifer` (reported to PopSeq). Those subprojects are not counted until it is fixed.
- Five `Population_Genomics` directories have no `config.txt` and are invisible to the inventory (PopSeq issue B4); `Scedosporium` holds several species and is not counted.
- The run table carries no `sample_alias` or `experiment_title`, which the 2026-10-26 survey used in `classify_study()`; the exclusion is therefore applied to strain, isolate and sample title only.
- Clonality is shown, not corrected: `largest_study_frac` flags dominance by one study, but nothing here clusters isolates by genotype (see `../sra_survey/pilot/clone_clusters.py` for that, per species).
