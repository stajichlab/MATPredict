# Cassette fields (`receptor_cassette_*`): before/after check

Baseline = frozen worktree run-pr-array-final (0bdde27, `receptor_array_*` only); candidate = frozen run-pr-array-cassette
(9780812). 37 distinct Basidiomycota genomes (the 33 of the regression panel plus the 6 Agaricomycete study-panel genomes, 2 overlap),
HPCC jobs 29533166 (base) and 29533167 (cand), `scripts/run_clade_panel.slurm`, list `panel39.tsv` (37 rows).
`compare_cassette.py` compares every report key including `receptor_array_*`; only `receptor_cassette_*` and `receptor_arrays_note`
may differ. `summary.txt`: 0 differences. `loci.tsv` (PR calls, with the loci columns), `arrays.tsv` (all arrays).
Full test suite: 1075 passed.
