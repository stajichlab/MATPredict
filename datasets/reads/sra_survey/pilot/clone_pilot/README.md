# Clone detection pilot (2026-10-27): can read sketches say which runs are the same strain?

Question: can a k=21 sourmash sketch of the 2-million-read subsample (scaled 500, k-mers seen at least twice) group the runs of one species into clones, so that MAT counts can be given per genotype instead of per run? Method and script: `../clone_clusters.py`; the same-strain ceiling: `../clone_calibration.py`; sketches written by the worker with `SKETCH=1` (`SKETCH_SPLIT=1` for the halves).

Data: 252 runs of six species (the 260 submitted minus 8 not at NCBI: 5 *P. roqueforti* from PRJEB20132 and 3 others), and 24 runs (4 per species) sketched twice from alternate halves of their own reads. Files here: `<species>.pairs.tsv` (identity estimate of every pair), `<species>.clusters.tsv` (clusters at 0.999), `calibration.tsv` (halves).

## Results
**Two halves of one run (same genome, each at half the run's coverage)**, 24 runs: identity estimate 0.9345 to 0.9993, median 0.9975. Only a minority reach 0.999, so a 0.999 cut-off calls most same-strain comparisons "different". One run (SRR4023060, 14,000 hashes, a low-yield set) gives 0.935.

**Between runs of one species vs the halves ceiling**

| Species | pairs | between median | p95 | max | halves min / median (4 runs) |
|---|---|---|---|---|---|
| *A. fumigatus* | 1,225 | 0.9956 | 0.9977 | 0.9990 | 0.9924 / 0.9989 |
| *C. posadasii* | 1,225 | 0.9941 | 0.9957 | 0.9989 | 0.9972 / 0.9984 |
| *F. oxysporum* | 1,225 | 0.9766 | 0.9892 | 0.9991 | 0.9758 / 0.9929 |
| *P. roqueforti* | 351 | 0.9920 | 0.9954 | 0.9988 | 0.9954 / 0.9955 |
| *T. reesei* | 630 | 0.9910 | 0.9994 | 0.9999 | 0.9345 / 0.9979 |
| *V. albo-atrum* | 10 | 0.8463 | 0.9987 | 0.9987 | 0.9979 / 0.9988 |

**Number of clones (single linkage) by identity threshold**

| Species (runs) | 0.999 | 0.998 | 0.997 | 0.995 | 0.990 |
|---|---|---|---|---|---|
| *A. fumigatus* (50) | 50 | 31 (largest 18) | 9 (largest 42) | 2 | 1 |
| *C. posadasii* (50) | 50 | 41 (6) | 41 (6) | 3 | 1 |
| *F. oxysporum* (50) | 49 | 48 (2) | 44 (4) | 43 | 27 |
| *P. roqueforti* (27) | 27 | 22 (6) | 22 (6) | 13 | 2 |
| *T. reesei* (36) | 23 (largest 14) | 21 (16) | 21 (16) | 20 | 3 |
| *V. albo-atrum* (5) | 5 | 4 | 4 | 3 | 3 |

## Reading
- *T. reesei* has a clear clone of 14 to 17 runs at every threshold from 0.999 to 0.995. That fits the many laboratory strains derived from QM6a. This is the case the method is for.
- *C. posadasii* and *P. roqueforti* show a plateau between 0.998 and 0.997 (41 and 22 clones, largest group 6): groups of near-identical runs, separated from the rest.
- *A. fumigatus* has no plateau: 31 clones at 0.998, 9 at 0.997 with one chain of 42 runs. The species is so uniform at the k-mer level that between-run similarity (median 0.9956) sits inside the same-strain range. Single linkage chains through it. Sketch similarity cannot separate clones from the rest here.
- *F. oxysporum* is diverse (median 0.977), so almost everything is a singleton at any strict cut-off; a handful of pairs (largest group 4) are near-identical.
- The ceiling depends on the run (coverage, yield): a single threshold for all species and runs is not defensible.

## Linkage test (2026-10-28)
The same pair tables clustered again (`clone_clusters.py --from-pairs`), as clones / size of the largest clone:

| Species (runs) | threshold | single | average | complete |
|---|---|---|---|---|
| *A. fumigatus* (50) | 0.998 | 31 / 18 | 44 / 3 | 45 / 2 |
| | 0.997 | 9 / 42 | 27 / 13 | 31 / 4 |
| | 0.995 | 2 / 46 | 5 / 43 | 15 / 11 |
| *C. posadasii* (50) | 0.998 | 41 / 6 | 42 / 5 | 43 / 4 |
| | 0.997 | 41 / 6 | 41 / 6 | 42 / 4 |
| *F. oxysporum* (50) | 0.998 | 48 / 2 | 48 / 2 | 48 / 2 |
| *P. roqueforti* (27) | 0.998 | 22 / 6 | 24 / 4 | 25 / 2 |
| *T. reesei* (36) | 0.998 | 21 / 16 | 22 / 15 | 24 / 11 |
| | 0.995 | 20 / 17 | 21 / 16 | 21 / 16 |
| *V. albo-atrum* (5) | 0.998 | 4 / 2 | 4 / 2 | 4 / 2 |

- **Chaining is the cause of the *A. fumigatus* failure, and average linkage removes it.** At 0.998 the 42-run chain becomes 44 clones, the largest of 3 runs (39 singletons, 4 pairs, one triple).
- **The other species do not depend on the linkage choice.** *C. posadasii* gives 41 to 43 clones at 0.998 and 0.997 for all three linkages; *T. reesei* keeps a clone of 11 to 16 runs; *F. oxysporum* has the same 48 clones.
- **Average linkage at 0.998 is the new default.** It does not stop merging inside the *T. reesei* clone at a high floor: its least similar pair is 0.996 (mean pair at least 0.998), so `min_ani` is reported per clone.
- **Still open: whether a 0.998 pair in *A. fumigatus* is a clone or two close strains.** Its between-run distribution has p95 0.9977 and p99 0.9983, so a few of the 1,225 pairs reach 0.998 by being close relatives, not copies. A per-run ceiling (a split sketch of every run) or an assembly-based identity of the same BioSamples would decide; neither is done.

## Decision this supports
Report **clone groups as a flag, not as a correction**: a species page may say "N of these runs fall in near-identical groups (identity at least 0.998)" only where the clone-count curve has a plateau, and must not give "per-clone" frequencies for species without one (*A. fumigatus* here). Use average linkage at 0.998 (done in `clone_clusters.py`). Not yet tested: a ceiling per run (needs a split sketch of every run), mash/skani on the assemblies of the same BioSamples as a cross-check, and whether more reads per run raise the ceiling.
