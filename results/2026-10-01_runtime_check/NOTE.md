# Runtime check: 076afe4 vs 52b3ff9 (2026-10-01)

## Question
The held-out rerun note reported a median of ~260 s/genome on 52b3ff9 against
~52 s on 076afe4. Is the slowdown real, and where does the time go? Accuracy
first: only call-neutral optimisations may be proposed.

## Data and code version
- Frozen worktrees .claude/worktrees/run-076afe4 and run-52b3ff9.
- Panel: 20 LCG Mucoromycotina genomes (14 called, 6 uncalled, random) + 5 Zygo
  genomes (panel.tsv). Command as the held-out runs:
  `detect --phylum Mucoromycota`, no taxid.
- One SLURM job (29323077), one node (r44), 16 CPUs, 8 genomes at a time,
  old/new order alternated per genome (tasks.tsv). Then cProfile on 3 genomes.

## Results
1. The 5x premise does not hold. The recorded wall_seconds of the two held-out
   runs give, for the same 621 LCG genomes, median 162 s (076afe4) vs 182 s
   (52b3ff9), ratio 1.24 (p10 0.84, p90 2.57); the runs were on different nodes.
2. Controlled, same node (timing_per_genome.tsv): median 55.1 s vs 64.9 s;
   per-genome ratio median 1.15 (range 0.96-1.40); total 1,406 s vs 1,622 s
   (+15%).
3. Where the time goes (cProfile, 3 genomes):

| Genome | Total old -> new (s) | tblastn search | exonerate polish | Python-side new steps |
|---|---|---|---|---|
| A. blakesleeana NRRL 1300 | 20.7 -> 25.7 | - | 8.9 -> 10.8 | ~0.8 |
| M. aromaticus NRRL A-17745 | 30.5 -> 38.8 | 14.9 -> 20.2 | 10.4 -> 13.4 | ~1.2 |
| Actinomucor NRRL A-23671 | 81.5 -> 95.4 | 7.6 -> 11.0 | 69.5 -> 79.3 | ~1.2 |

   - Over 90% of wall time is in external tblastn and exonerate subprocesses in
     both versions. Exonerate call counts are unchanged (10 and 18); each call
     aligns more reference proteins.
   - Cause: the Mucoromycota reference set grew from 41 to 51 proteins
     (21,767 -> 26,590 aa, +22%) with the Umbelopsis Plus/Minus and
     S. racemosum NRRL 2496 records (merged a2fe1b4). tblastn and exonerate
     time scale with query size, which matches the ~15-24% rise.
   - New Python steps are small: classifier scoring 0.3-0.5 s, HMM loading
     0.2-0.3 s, MAT-gene gate 0.1 s, V3 ranking 0.5-0.8 s, two_idiomorphs and
     split-locus ~0. CAAX scan does not run for Mucoromycota (0.0 s).

## What changed in detection
Nothing; measurement only. No code committed.

## Call-neutral optimisation candidates (not prototyped)
- Cache parsed YAML records/roster (safe_load ~2.5 s/genome, 124-129 loads):
  ~4-7% per genome; inputs identical, so call-neutral by construction.
- Load classifier HMMs once per process (0.2-0.3 s/genome).
- Not proposed: reducing tblastn/exonerate query sets or adding tblastn
  threads; those change search inputs or hit order and are not provably
  call-neutral.
The slowdown is the cost of a larger, more accurate reference set; no change
to it is recommended.

## Limits
25 genomes, one node, 8-way concurrency; profiles on 3 genomes.

## Files
panel.tsv, tasks.tsv, run_timing.slurm, time/timings.tsv,
timing_per_genome.tsv, prof/*/<genome>/profile.pstats, logs/.
