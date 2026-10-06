# Receptor arrays (`receptor_array_support`): before/after check

Branch `pr-array-support`; baseline = `origin/main` 2d1860a, candidate = this branch. Report-only change.

- `compare_calls.py`: reads two sets of `detection_report.yaml`, checks every `detected` entry and every top-level key
  (bar the two new ones) is identical once the `receptor_array_*` fields are dropped, writes `loci.tsv` with the five array
  columns (`receptor_arrays.loci_columns`) and a summary.
- `summary_final37.txt`, `loci_final37.tsv`: final code (df08f63) on the 37 Basidiomycota genomes (the 33 of the regression
  panel plus 6 Agaricomycete study-panel genomes, 2 overlap). 57 calls on each side, 0 differences outside the array fields.
- `regression_summary.md`: `scripts/run_regression_panel.sh` (163 genomes + Zygo 23 on both inputs), earlier candidate commit
  1d14348 (arrays without the 50%-coverage locus filter): 0 loci changed in every panel.
  `summary_regression_oldcand.txt`: the same reports through `compare_calls.py`, 0 differences.
- `dump_hits.py`, `dump_hits.sh`: dumped the pre-clustering PR hits of two genomes; Trametes versicolor has 372 receptor
  HSPs, 67 merged loci, 6 with >= 50% reference coverage, which set the locus filter.

## Final run (rename to `receptor_array_*`, commit 0bdde27)

Frozen worktrees run-pr-array-base (main 2d1860a) vs run-pr-array-final (0bdde27); `scripts/run_regression_panel.sh`
on the 163-genome panel plus Zygo 23 (HPCC `results/2026-10-06_regression_pr_array_final/`). Loci changed: 0 in
ascomycota (30), basidiomycota (33), mucoromycota (80), record_sources (23) and zygo23 (46 reports). `compare_calls.py`
per panel: 0 differences outside the `receptor_array_*` fields; PR calls with array fields: 14 (basidiomycota), 2
(record_sources). Summaries in `final_panel/`. Full tests: 1065 passed (mafft 7.505).
