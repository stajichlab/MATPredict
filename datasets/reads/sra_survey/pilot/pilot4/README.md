# Pilot test 2: 80 runs, four species, checked against `detect` (2026-10-26)

20 runs each of *Pyricularia oryzae* (no reference in the database), *Botrytis cinerea*, *Verticillium dahliae* and *Aspergillus flavus* (genus records only), 8 at a time on the HPCC (`../pilot4_runs.tsv` is the run list; `SUMMARY.md`, `summary_by_group.tsv`, `results_all.tsv` and `failures.tsv` are the results; final: 79 of 80 runs done). Each panel is the MAT1-1 and MAT1-2 locus with flanks cut from two BFD assemblies on the `detect` coordinates, as for *N. crassa* (`panels` on the cluster: `results/2026-10-26_v070_submit/panels/`; to be committed with the species' panel folders).

**Truth.** The `detect` idiomorph of a BFD assembly that shares the run's BioSample (via NCBI's assembly summary). 53 of the 79 finished runs have one.

## Result
- **All 51 typed runs with a truth agree with `detect`: 16/16 *A. flavus*, 11/11 *B. cinerea*, 14/14 *P. oryzae*, 10/10 *V. dahliae*.** Two *P. oryzae* runs with a truth were not called (`low_depth`, below). No run was called with the wrong idiomorph.
- Typed calls without a truth are plausible: *V. dahliae* 8 MAT1-2 and 0 MAT1-1 among the runs without truth (the 10 with truth: 9 MAT1-2, 1 MAT1-1, all right) (BFD assemblies: 54 MAT1-2, 5 MAT1-1); *A. flavus* 3 MAT1-2 and 1 MAT1-1; *B. cinerea* 5 MAT1-2, 2 MAT1-1.
- The other idiomorph appears at low depth in many runs (the `trace_` flag). Among typed runs the other idiomorph's depth is a median 0.000, 90th percentile 0.19 and at most 0.41 of the shared-flank depth.

## Calls that are not a simple one-idiomorph call
| Run | Species | Call | What it looks like |
|---|---|---|---|
| SRR23388381 (UFVPY171), SRR29065219 (UFVPY657) | *P. oryzae* | `low_depth`, truth MAT1-1 and MAT1-2 | Shared-flank depth 0.3 and 0.26 although the runs hold 6 and 13 million spots: the panel flank barely matches. A different lineage or host-specialised group from the rice-blast strains the panel was cut from is the likely cause (not checked). A second panel would be needed. The pipeline declined to call them. |
| SRR14970071, SRR14970059 (VD991) | *V. dahliae* | `none` | Neither idiomorph has any depth (0.00) while the shared flank is at 2.2x. Not explained: a divergent MAT locus, a different lineage or species, or a deletion. No truth. |
| SRR32366559 ("Bc116 x Bc448") | *B. cinerea* | `both` | MAT1-1 5.9x and MAT1-2 2.3x, both breadths high (1.00, 0.82). The name says a cross: a mixed or heterokaryotic sample is plausible. |
| SRR24910575 (Bcvcg1) | *B. cinerea* | `both` | MAT1-2 4.6x, MAT1-1 0.9x (ratio 0.2, flagged), breadths 0.89 and 0.58. Borderline between a trace and a real second idiomorph. |

## Not completed
- `ERR2188722` (*P. oryzae*): **unavailable at NCBI.** Not in NCBI's AWS copy, and NCBI itself answers "no data (404)" for it (`prefetch`, `srapath`), while ENA has the FASTQ. About 1 run in 100 so far (the two other ERR runs, in the *N. crassa* test, were in the AWS copy). The worker now records `unavailable_at_ncbi`. Fetching such runs from ENA (0.13 to 0.9 MB/s) is possible but slow; not done.
- `SRR8238192` (*V. dahliae*): a **128 Gb run (429 million spots), 39.5 GB to download**. It finished (download 36 minutes, conversion 17 s, typing 63 s) and was typed MAT1-2, agreeing with `detect`. It was chosen by my rule (the largest paired run per sample), which was wrong: the next selection takes the run nearest a few Gb and caps the size. The windowed conversion made the size irrelevant for time after the download.

## Speed (measured)
175 GB for the 79 finished runs (39.5 GB of it one run); median 158 s per run; 12 to 15 MB/s per stream with 8 streams at once (about 110 MB/s together), against 17 to 20 MB/s per stream at 4 at once: it scales, with some loss per stream.

## Next
- Look at the four non-simple calls (a second *P. oryzae* panel from another lineage; the VD991 locus; the two *B. cinerea* samples).
- Choose runs by size (nearest 3 Gb, at most about 10 Gb).
- Panels and the same test for *Z. tritici*, *S. sclerotiorum* (homothallic control) and *C. parasitica*.
