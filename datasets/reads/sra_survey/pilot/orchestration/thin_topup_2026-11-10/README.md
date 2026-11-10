# Top-up of the thin-evidence species (prepared 2026-11-10, to run after the HPCC shutdown)

PR #327 typed 17 species that used to be held for thin evidence (one assembled strain, or 1 to 3 runs); they are UNVALIDATED with `single_locus` or `few_runs_<n>`. They were typed only from their plan's run list (at most 50 runs). This folder tops them up to the cap of 700 runs that the other species have (round 3).

- `species.txt`: the 17 species, with eligible and typed counts from the 2026-11-06 snapshot. About 340 runs are expected, mostly *F. sambucinum*, *P. fijiensis*, *Madurella*, *C. amylolentus* and *F. andiyazi*; the `few_runs` species have 4 to 9 eligible isolates, so they gain little unless PopSeq's newer table has more.
- `run.sh refresh | lists | submit | check`: new availability snapshot, run lists from `expand_runs.py --species-file species.txt --cap 700 --per-study 700`, stage `thin_topup` per wave, then a list of runs to rerun.
- Merge: `append_stage.py` per wave. It appends the new runs, leaves existing rows alone, and decides the topped-up species again. The thin-evidence cap stays: more runs do not lift `single_locus`, because the panel still comes from one strain.

Choices to confirm: no per-study limit (as in `expand300b`), and the 700 cap. The snapshot the refresh writes is committed with the results.
