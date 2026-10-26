# Pilot test 1: 20 *Neurospora crassa* runs (2026-10-26)

`../stream_type_sra.slurm` on the HPCC (`exfab`, 4 tasks at a time, 4 CPUs and 8 GB each) against the panel in `panels/neurospora_crassa/`. `../ncrassa20_runs.tsv` is the run list (ENA accessions and sizes); `results.tsv` is one row per run. Reads are fetched as `.sra` from NCBI's copy on AWS open data, cut into 8 evenly spaced windows of spots (2 million spots in all), typed, and deleted; only the small result files stay.

## Result
- **Every paired run typed to one idiomorph (16 of 16): 9 MAT1-1 and 7 MAT1-2.** The two strains whose mating type is in the strain name were both right: OR74A (FGSC2489, mat A) typed MAT1-1 and 74-ORS-6a (FGSC4200, mat a) typed MAT1-2. The Louisiana wild isolates split 4 MAT1-1 and 2 MAT1-2.
- **The four single-end runs from the "*Neurospora crassa* Variation" study (PRJNA355077) were not called**: two `none` (idiomorph depth 5.4 and 7.7 but breadth only 0.04 to 0.06, where 5x depth would give about 0.99) and two `low_depth` (shared flank depth 0.08 and 0.12). The depth is concentrated on a few positions, which is the pattern of reduced-representation data (RAD or GBS-style libraries) filed as WGS, not of a whole genome. This is inferred from the depth and breadth only; the read lengths and library were not inspected. The pipeline did not call them wrongly. Single-end WGS runs are not shown to work either way by this test.
- **Every typed run has a `trace_` flag** for the other idiomorph, almost all at 0.0 to 0.5x. One is higher: SRR28962352 (isolate 7833) has MAT1-2 depth 3.99 against MAT1-1 depth 7.83 (breadth 0.26). Not explained: panel cross-talk, a mixed sample or a heterokaryon are all possible.

## Speed and transfer (measured)
- 22.7 GB of `.sra` for 20 runs (median 0.9 GB). Download 17 MB/s on average per stream (3 to 227 s per run); conversion 11 to 39 s; typing 14 to 79 s. Median task 108 s.
- The 19-run array finished in 10.7 minutes of wall time at 4 tasks at a time. Peak memory 4.9 GB.
- Same cluster, same day: ENA's FTP gave 0.13 MB/s on one stream; NCBI's `sra-instant` path returned 404; the AWS copy gave 20 to 26 MB/s per stream and 76 MB/s with four.
- The two European runs (ERR11702005, ERR11702008) are in the AWS bucket too.

## Next
- Look at the 7833 trace and at what the single-end "Variation" reads are.
- Same test on the other pilot species once their panels exist (*P. oryzae*, with *C. parasitica* next, then *Z. tritici*, *B. cinerea*, *V. dahliae*, *A. flavus*, *S. sclerotiorum*).
