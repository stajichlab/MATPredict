# tblastn `-seg` policy for the slowest genomes (#276)

Status: opt-in, default unchanged. Written 2026-11-06. Measurements are in the #276 thread; files are under
`results/2026-11-06_*` on the HPCC (not committed).

## Problem
The genome-wide localization runs one batched `tblastn -seg no` with every reference protein of the routed families
(`search_localize`). For `tblastn` the BLAST+ default is `-seg 12 2.2 2.5` (low-complexity filtering ON), so production does less
filtering than stock. On most genomes this costs seconds, but on a few it costs more than an hour: in the v0.8.1 run 13
Ascomycota:MAT genomes took over 600 s and two (*P. ulei* IBUN, DampeI2T2) hit the 3,600 s limit and needed a rerun.

## What is slow
Per reference protein, one CPU, `-seg no`, on an EPYC node: SLA2 (17 proteins) 1,090 s on DampeI2T2 and over 1,500 s on IBUN;
MAT1-1-1 (22) 386 s and 687 s; APN2 (19) 184 s and 334 s; MAT1-2-1 (26) 7 s and 11 s; the other nine gene groups 0 to 4 s. A
median-speed genome takes 6 to 11 s for each of those three groups. The cost is specific to particular proteins: the fraction of
a query that `segmasker` marks is 6 to 17% for the slow queries and up to 40% for fast ones.

## What does not work
- Plain `-seg yes` on every query: 188 s to 1 s on the slowest query, but the bitscore of the real hit changes (1424 to 1175) and
  calls change. On a 300-genome Ascomycota panel 231 genomes were identical, 7 changed family, idiomorph or contig; on a
  120-genome Basidiomycota panel it loses the short pheromone precursors (MFalpha 23.9 bits, pheromone_B43 25 to 27 bits,
  fungal_mating_type_pheromone 25 bits) and 2 calls change.
- It gives no speed-up on typical genomes: same node (r12), 60 Ascomycota genomes, 16 at a time: median 82 s with `-seg no`,
  86 s with soft-masked SEG (1.05x). Earlier ratios of 0.47x compared against baseline runs on older nodes and are not valid.

## What works
`-seg yes -soft_masking true` (masked regions used for seeds only; alignments scored and extended through them) on the long query
only. Two ways to choose the queries, both opt-in through environment variables and both implemented as two `tblastn` calls on the
same database, with the output merged in the reference file's query order:

- `MATPREDICT_TBLASTN_SEG_GENES=SLA2,APN2,MAT1-1-1` (recommended): only these gene names get the filtered search.
- `MATPREDICT_TBLASTN_SEG_MIN_AA=250` (default when only `MATPREDICT_TBLASTN_SEG=yes` is set): queries of at least that many residues.

Results against the v0.8.1 reports (frozen code, `--preset v0.8.1`, EPYC nodes):

| Policy | Ascomycota 300 identical | Basidiomycota 120 identical | Mucoromycotina 90 identical | Outliers IBUN / DampeI2T2 |
|---|---|---|---|---|
| `-seg no` (production) | 300 | 120 | 90 | 3,063 s / 1,685 s |
| soft-masked SEG on every query | 282 (0 call changes) | 96 (2 call changes) | 89 | not run |
| length rule, 250 aa | 290 (0 call changes) | 105 (0 call changes) | 90 | 468 s / 271 s, loci identical |
| gene list | 296 (0 changes of family, idiomorph, confidence, gene set) | 120 | not run | 468 s / 288 s, loci identical |

With the gene list, the four remaining Ascomycota differences are locus spans: two of 9 and 27 bp, and two of 12 to 20 kb at a locus
edge with an identical gene set and identical scores.

## Not established
- Only the Ascomycota MAT family was timed per gene. Mucoromycotina and the yeast families have no slow tail in the v0.8.1 run
  (maximum 63 to 126 s). Basidiomycota has one (HD up to 3,093 s, PR up to 1,891 s); its slow proteins are being timed.
- Hardware matters more than this setting for the whole run: the v0.8.1 chunks on the older AMD Opteron nodes ran about 2.15 times
  slower than on EPYC nodes (use `--constraint=milan` for timing).
