# tblastn-fast as the localize engine: speed, sensitivity and call changes
Status: open (curator review of the two lost calls and the core-model changes pending)

## Question
Does `tblastn -task tblastn-fast` change speed, hit sensitivity or final calls
in the localize step of `detect`, compared with the default `tblastn`?

## Data and code version
- Code: `73f56b8` (base, default task, frozen worktree `.claude/worktrees/run-73f56b8`)
  and `9b15d43` (candidate, only `DEFAULT_LOCALIZE_TASK = "tblastn-fast"` differs,
  frozen worktree `.claude/worktrees/run-9b15d43`). Both on branch `localize-engine`
  (local; `9b15d43` is an EXPERIMENT commit and must not be merged as is).
- Database: same tree as the code. Both commits have the same `db/`.
- BLAST: tblastn 2.17.0+ in both the benchmark (`module ncbi-blast/2.17.0+`) and the
  panel jobs (pixi default env).
- Inputs: (a) benchmark: 447 curated reference proteins (`build_reference_fasta`
  over the whole `db/`, count from `analysis/tblastn_fast/run.log`; the "181" in
  the `reference_fasta.py` docstring is out of date), 15 truth
  assemblies from `results/truth_assemblies.json`; (b) regression panel:
  `testset/regression_panel.tsv` (163 genomes) plus Zygo 23.

## Method
1. Benchmark: for each of 15 assemblies, run `tblastn` and `tblastn -task tblastn-fast`
   with `-seg no -num_threads 4`, same query, same database, 2 replicates each
   (alternating). Script: `analysis/tblastn_fast/run_bench.sh`.
2. Compare HSP counts, time, truth placement overlap and best-hit retention
   (`analysis/tblastn_fast/analyze.py`).
3. Characterize default best hits (e <= 1e-5) that fast mode does not reproduce
   (`analysis/tblastn_fast/lost_hits.py`, table `lost_hits.tsv`).
4. Regression panel: `scripts/run_regression_panel.sh` on the two frozen worktrees,
   partition `short`, 16 CPUs per job, jobs 29893802-29893814. Diff by
   `scripts/regression_check.py`.
5. Trace the two lost calls in the per-genome `detection_report.yaml` files.

## Results
### Benchmark (15 assemblies, median of 2 replicates)
| | default | tblastn-fast |
|---|---|---|
| Total time | 775 s | 164 s (4.7x; per genome 3.1-5.4x) |
| HSPs at e <= 10 | 36,699 | 14,031 |
| Truth gene placements overlapped | 41/44 | 41/44 |
| Default best hits (e <= 1e-5) kept at the same place | | 2,881/3,326 (86.6%) |

Source: `analysis/tblastn_fast/out/timing.tsv`, `analyze.py` output. The 3 truth
placements missed are all on GCF_000182925.2 (*N. crassa*) in both modes. The cause
is not investigated.

### Hits lost by fast mode (445 default best hits, e <= 1e-5)
- Identity: median 27.6%, maximum 48.8%. Query coverage median 43%.
- Role: 363 core_MAT, 66 flanking_conserved, 16 flanking_variable.
- Bitscore: all >= 40.8 (median 59.7); 109 have >= 80 bits; 222 have e <= 1e-10.
- 338/445 (76%) of the loci are still found by another query of the same gene name.
- 107 are not recovered by any same-gene query: 62 core_MAT, 29 flanking_conserved,
  16 flanking_variable. Largest groups: rnhA 26, glrA 12, STE3v2 9, STE3a2 8.
- All 15 genomes are Ascomycota. No Mucoromycota or Basidiomycota genome was in
  the benchmark.

### Regression panel (`results/2026-10-12_regression_tblastn_fast/diff/`)
`diff/summary.md` lists 6,592 changed loci. `regression_diff.tsv` has 7,648 rows
over the five panels (the two counts use different units). Of those rows, 5,991
carry `withheld_reason_changed`: mostly withheld candidate clusters that vanish or
change reason. Rows that touch a called locus on either side: 206.

| change type, call-touching rows | count |
|---|---|
| call_lost | 2 |
| call_gained | 0 |
| idiomorph_changed | 0 |
| confidence_changed | 6 (4 high->medium, 2 medium->high) |
| locus_class_changed | 4 |
| classifier_shift (>= 5 bits) | 3 (Zygo 23, Minus *Cunninghamella*, margin 46-50 -> 34.5, call unchanged) |
| core_model_changed | 21 |
| gene_set_changed / span_changed | 24 / 48 |

Zygo 23: no call lost or gained, no idiomorph change.

### The two lost calls (both Basidiomycota, both `undetermined` idiomorph, medium confidence)
1. **GCA_023273805.1, family Basidiomycota:MAT, CP096880.1:1,935,834-1,981,717.**
   Base and candidate both find FAO1 and RPL39. The candidate scores
   `fraction_found` 0.333 against the 0.5 floor and withholds it
   (`below_fraction_floor`). The base reported it with
   `ambiguous_with: [Basidiomycota:Balpha, Basidiomycota:bLocus]`. The pipeline
   reports a below-floor cluster when two or more families clear the floor
   (`is_ambiguous`, `pipeline.py` around line 2915). In the candidate the Balpha
   hit (bap3-3) at this locus is not reported and the bLocus cluster is shorter
   (to 1,964,900), so the rescue does not apply.
2. **GCA_056247695.1, family Basidiomycota:PR, JBVYWQ010001089.1:920-8,651.**
   The base found three genes: pheromone_receptor (exonerate, 75.9%, 248 bits),
   a CAAX scan precursor, and `fungal_mating_type_pheromone` by tblastn at 37.0%
   identity, 26.9 bits, e = 0.62. The candidate does not report the tblastn
   pheromone hit. `fraction_found` falls to 0.25 and the locus is withheld. The
   base call was already `unverified` (CAAX scan, `idiomorph_gene_only`).
   The reference pheromone queries are 25-43 aa long.

In both cases the base call rested on evidence at or below noise level. The
data here do not show whether either base call was a true MAT locus.

## What changed in detection
Branch `localize-engine`, commit `73f56b8`: `search_localize(..., task=None)`
accepts `tblastn` or `tblastn-fast`; the default adds no flag. Tests:
`tests/detect/test_search.py` (26 pass). Commit `9b15d43` flips the default and is
for this experiment only.

## Limits
- One run of each side on a shared cluster. Job wall times (candidate 9:58-15:19,
  base 14:05-21:28 per group) include startup and queue contention. The 4.7x figure
  is the search step alone.
- The panel has 163 genomes and 46 Zygo 23 call rows. Two lost calls is a count on
  this panel, not a rate for the BFD/NCBI cohort.
- No ground truth was used for the two lost calls or for the 107 unrecovered hits.
- mmseqs is not measured yet.
- The Opus review `engine_scoring_implications.md` is model output. Its findings
  are not all re-verified. Treat its line references as leads.
- Short precursor queries (25-43 aa) may be the weakest case for the larger word
  size. This is a hypothesis from case 2. It is not measured.

## Curator decisions
None yet. Questions for the curator are in `open-questions.md`.

## Files
- `analysis/tblastn_fast/` (`run_bench.sh`, `analyze.py`, `lost_hits.py`,
  `lost_hits.tsv`, `out/timing.tsv`, `out/*.tsv.zst`)
- `analysis/engine_scoring_implications.md`
- `results/2026-10-12_regression_tblastn_fast/diff/summary.md` and per-panel
  `regression_summary.md`, `regression_diff.tsv`
