# Cassette fields (`receptor_cassette_*`): before/after check

Baseline = frozen worktree run-pr-array-final (0bdde27, `receptor_array_*` only); candidate = frozen run-pr-array-cassette
(9780812). 37 distinct Basidiomycota genomes (the 33 of the regression panel plus the 6 Agaricomycete study-panel genomes, 2 overlap),
HPCC jobs 29533166 (base) and 29533167 (cand), `scripts/run_clade_panel.slurm`, list `panel39.tsv` (37 rows).
`compare_cassette.py` compares every report key including `receptor_array_*`; only `receptor_cassette_*` and `receptor_arrays_note`
may differ. `summary.txt`: 0 differences. `loci.tsv` (PR calls, with the loci columns), `arrays.tsv` (all arrays).
Full test suite: 1075 passed.

## Rename (follow-up)

`receptor_cassette_caax_orfs` was renamed `receptor_cassette_max_caax_orfs`: it is the maximum number of strict-CAAX ORFs within the
window of any single locus in the array, not a count of cassettes (it can be 1 on an array with no cassette). The run outputs above
were produced under the old name; only the column header of `arrays.tsv`/`loci.tsv` and `compare_cassette.py` were relabeled, the
values are unchanged. The before/after was rerun on the final commit (see below).
