# Polish-cap rank: the mixed rank does not pass, and why the cap loses real loci

2026-09-26. Curator ruling: test a mixed rank (top 5 clusters by gene count +
the single highest-identity cluster) by replay; adopt only if it loses
nothing that either single rank keeps. **It does not pass. The rank was not
changed.**

## Method

`rank_simulation_mixed.py` replays a per-family cap of 6 on UNCAPPED runs
(evidence-diagnostics `evidence` rows). A call is lost when no kept cluster
of its family overlaps it. Panels: the 561-genome cap panel (`off/`), the
early-diverging Mucoromycota / Mortierellomycota / Kickxellomycota scans,
the 2,368-genome Serinales scan, and a new uncapped run of the four
Dothideomycetes genomes the round-2 curation lost
(`mixed_rank_dothideo/`, code+db `run-52cd292`, `--max-polished-clusters-per-family 0`).

## The first replay was wrong, and that found the cause

The plain replay (gene_count from the evidence row) said genes-first kept
all three Dothideomycetes calls, but the real capped run lost all three.
The real cap ranks on LIVE distinct genes (`_polish_rank` excludes
superseded hits); the evidence row counts them. A real MAT locus usually
draws a cross-hit from the other idiomorph's reference, which the first
resolution supersedes. Example, Zymoseptoria brevis GCA_000966595.1: the
true locus (MAT1-1-1 95.6%, with a superseded MAT1-2-1 cross-hit at 44.1%)
has 3 genes in the row but 2 live, so it ranked below six 3-gene clusters at
32.7-42.9% identity and was never polished.

`--live` subtracts each contig's pre-polish resolution losers. It reproduces
the real Dothideomycetes losses exactly (3/3) and gives 5 on the 561 panel
(real cap6 run: 4; the replay is approximate where a contig holds several
clusters).

## Results, live replay, N = 6 (`rank_simulation_mixed_live.txt`)

| panel | calls | genes-first (current) | identity-first | mixed (5+1) | mixed 4+2 |
|---|---:|---:|---:|---:|---:|
| Dothideomycetes 4 | 3 | 3 | 0 | 0 | 0 |
| cap panel 561 | 613 | 5 | 0 | 3 | 3 |
| Mucoromycota 293 | 242 | 0 | 2 | 0 | 0 |
| Mortierellomycota 100 | 6 | 0 | 0 | 0 | 0 |
| Kickxellomycota 190 | 7 | 1 | 0 | 1 | 1 |
| Serinales 2,368 | 2,647 | 5 | 3 | 5 | 3 |

Mixed loses 3 (cap panel), 1 (Kickxellomycota) and 2 (Serinales) calls that
identity-first keeps. Identity-first loses 2 Mucoromycota Minus calls that
genes-first keeps. No rank tested loses nothing the others keep.

## A candidate the curator may want to rule on (not implemented)

Ranking on ALL distinct genes, superseded included (the plain replay), loses
7 calls across these panels against 14 for the current live rank. It loses
only one call the live rank keeps (GCA_029290875.1, Ascomycota:MAT) and keeps
13 the live rank loses, among them all three Dothideomycetes cases, Bauco1,
Acidiella bohemica and Zymoseptoria brevis. It fails the strict "loses
nothing" bar by one call. Whether a superseded cross-hit should count toward
a cluster's rank is a curator decision.
