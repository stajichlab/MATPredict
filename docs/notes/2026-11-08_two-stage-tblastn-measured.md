# Two-stage `tblastn` search, measured: not faster, not identical (#287 step 2)

Date: 2026-11-08. Status: **not merged.** The code is on the branch `two-stage-search` (opt-in setting `MATPREDICT_TBLASTN_TWOSTAGE=yes`, tests
included) and stays there. Data: `results/2026-11-08_twostage*` on the HPCC.

## The idea

`tblastn` takes 62-64% of the wall time of a typical Pezizomycotina genome (#287 profile). The flank genes (SLA2, APN2, COX13) are about 45% of the query
residues. Search the core genes genome-wide, then the flank genes only in windows around the core hits, with the genome's full size given as `-dbsize`
so E-values match; a genome with no core hit gets the flank genes searched genome-wide. Estimated saving: up to about 40% of wall time.

## What was built

`search_localize` gains an opt-in two-stage branch (stage 1: `core_MAT` queries genome-wide; stage 2: the other queries against a database of
windows around every hit found, grown until the flank hits stop changing, 60 kb padding above the largest configured cluster gap of 50 kb).
It falls back to the plain search when a side is empty, the genome size cannot be read, or the windows do not reach a fixed point. 15 tests, full suite
1,904 passed.

## Result 1: it is slower

Same node, one allocation, 40 stratified Ascomycota genomes, four passes in the order plain, two-stage, two-stage, plain (so node type and warm-cache
effects cancel). Node c12 (AMD Opteron 6376; an older node, so absolute times are long, but every pass used the same one).

| | total (s) | median per genome (s) |
|---|---|---|
| plain | 3,591 and 3,300 | 94.5 and 82.0 |
| two-stage | 4,792 and 4,500 | 80.5 and 75.0 |

Per genome, two-stage / plain (ABBA mean): **median 1.17, total 1.35**; faster for 9 of 40 genomes, fastest 0.76 (24% saved), slowest 3.12.
(The first runs of the panels, on the `short` partition against milan baselines, gave ratios of 1.5-1.8 and were invalid for timing: different node types.)

Likely causes (not separately measured here): `tblastn` translates the whole subject on every call, so each extra call adds a fixed cost that caps what
splitting the query set can save; a round that confirms the windows did not change always follows the last growth round; and the many spurious core
hits (HMG-box paralogs, 15-35 withheld clusters per genome) each make a window.

## Result 2: it is not identical

Against the v0.8.1 reports (frozen database, `--preset v0.8.1`), calls (family, contig, coordinates, idiomorph, tier, genes):

| Panel | Genomes | Calls identical | Withheld-locus lists identical |
|---|---|---|---|
| Ascomycota | 300 | 300 | 286 |
| Basidiomycota | 120 | **117** (3 differ) | 87 |
| Mucoromycotina | 90 | 90 | 83 |

- A fixed 100 kb window (first version) lost the far end of a 270 kb *Cryptococcus*-type MAT locus (GCA_006352305.1: locus start 538,707 became 624,794).
  Growing the windows did not fix the Basidiomycota genomes: three still differ; why is not determined here.
- The withheld lists lose flank-only clusters on contigs that have no core hit (for example APN2, COX13 and SLA2 alone). They are withheld, never calls,
  but no windowed search can recover them without a genome-wide flank search, which removes the saving.

## Conclusion

The two-stage search fails both bars set for #287 (no change in calls; a speed-up on a same-node comparison). It is not worth more work as designed.
The levers that remain are those already measured: the choice of node (`--constraint=milan`; the older nodes are about 2.15 times slower) and the
YAML-loading change (#292, about 7% of a run). A different approach would be to avoid the repeated whole-genome translation, for example by giving
`tblastn` a pre-translated database; that is a separate experiment and was not tried.
