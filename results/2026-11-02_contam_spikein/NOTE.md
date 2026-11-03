# How much contamination can GmG see in a read set? (spike-in, 2026-11-02)

Method and scripts: `datasets/reads/contamination_spikein/`. Real read sets mixed at known fractions, screened with GmG 0.0.1
(index 2026.10.5, `--max-reads N` = N pairs). Two backgrounds (*C. posadasii* GT-134, *C. immitis* 3M3), each mixed with
*Histoplasma* (NM_7898) and *Lecanicillium* (CA10) at 0.5, 1, 2, 5, 10 and 25% of the pairs, and with the other *Coccidioides*
species as a same-genus negative. 38 mixes, four depths: 152 screenings. Tables here: `spikein.rows.tsv` (every screening),
`spikein.detection.tsv` (per cell), `spikein.false_positives.json`, `flagged_in_clean_552.tsv`.

## Detection: GmG's own `mixture_suspected` flag

Smallest depth at which the flag was raised for both backgrounds (2 of 2); n = 2 per cell, so treat the edges as approximate.

| contaminant pairs | *Histoplasma* | *Lecanicillium* |
|---|---|---|
| 0.5% | 6M | 6M |
| 1% | 2M | 6M |
| 2% | 2M | 2M |
| 5% | 500k | 500k |
| 10% | 500k | 500k |
| 25% | 100k | 100k |

So the 2M pairs the reads pipeline streams reliably show a contaminant at about 2% of the pairs and sometimes 1%; 6M reaches 0.5%.
Estimated genome coverage of a pure control at 100k / 500k / 2M / 6M pairs: 0.6-0.7x / 2.2-3.2x / 8.9-13x / 31-44x.

## Specificity

- 8 pure controls (two backgrounds, four depths): no flag, at any depth.
- 48 same-genus mixtures (*C. immitis* and *C. posadasii*): no flag. In one (*C. immitis* with 25% *C. posadasii* pairs, 6M) GmG's species call changed to *C. posadasii*. A genus-level screen does not see this kind of mix, as expected.
- Every set shows a "second genus" (*Aspergillus*, *Fusarium*, *Colletotrichum*) from chance k-mer matches; it has 3-12 hashes and is not flagged.

## A rule on hashes and concentration did not do better

A candidate rule (second genus named, at least 5 hashes, at least 70% of them in one reference genome) fires earlier in some cells
(5% *Histoplasma* at 100k pairs) but not reliably: it misses 5% *Histoplasma* at 2M and 2% at 6M, because the contaminant's hashes
spread over several genomes of the genus as depth grows. GmG's own flag was monotonic in both depth and fraction here. Use GmG's flag.

## The 552 "clean" *Coccidioides* sets are not all clean

The TE project's species audit passed 552 of 558 sets (1 mixture, 3 wrong-genus, 2 wrong-species excluded). Screened at 2M pairs,
GmG's flag fires on 19 of them (`flagged_in_clean_552.tsv`); the primary genus is *Coccidioides* in all 19:

| second genus | sets | share of hashes it explains | note |
|---|---|---|---|
| *Histoplasma* | 4 (M192, NM_0071, M158, M216) | 31-37% (1,569-1,985 hashes) | its hashes cover 39-57% of a *Histoplasma* genome: heavy contamination |
| *Candidozyma* | 7 (one at 20%: B12526) | 1.1-20% | the six others are at 1.1-2.9% |
| *Calonectria* | 1 (M223) | 7% | |
| *Malassezia* | 3 | 2.2-2.6% | |
| *Aspergillus* | 2 | 2.4-2.6% | |
| *Nakaseomyces* | 2 | 1.2-1.5% | |

The thirteen at 1-3% are at the detection limit measured above and are not confirmed: read-level assignment (mapping the reads to the
named genome) should confirm them. The four *Histoplasma* sets and B12526 are far above the limit.

## Limits

n = 2 per cell; one pair of backgrounds and two contaminants; the mix takes each source's reads in file order; read length differs
between sets (93-150 bp), so equal pairs are not equal bases; GmG 0.0.1 and one index release.
