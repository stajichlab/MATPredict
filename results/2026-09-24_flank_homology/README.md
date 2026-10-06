# SLA2 / APN2 homology in Pezizales, from annotated proteomes

Question (curator, 2026-09-24): can SLA2/APN2 from relatively close species
build a homology case in Pezizales genomes, to support using them as MATtub
flanks?

Inputs: the 53 of 172 BFD Pezizales assemblies that have an NCBI protein
annotation (`annotated.txt`, `proteomes/`). Queries (`queries.faa`): Morchella
importuna SLA2 (end4, AVI60811.1) and APN2 (AVI60802.1+AVI60803.1, annotated as
two CDS) from the locus deposit KY782629.1, plus A. fumigatus SLA2/APN2 from
the curated record as a distant Pezizomycotina comparison. MAT1-1-1/MAT1-2-1
from both sources locate the MAT gene. Reciprocal check: each best hit blasted
back against the A. fumigatus Af293 proteome, where SLA2 = XP_754991.1 and
APN2 = XP_754988.2. Positions from each assembly's GFF. Scripts: `run.py`,
`summarize.py`; per-genome table in `summary.txt`, raw in `besthits.tsv`.

## Homology

| gene | found | reciprocal best hit to Af293 | E, Morchella query | E, A. fumigatus query | second protein hit |
|---|---:|---:|---|---|---|
| SLA2 | 53/53 | 53/53 | median 0, worst 1e-57 | median 0, worst 1e-85 | none in 49 genomes |
| APN2 | 46/53 | 46/46 | median 2e-170, worst 3e-91 | median 1e-166, worst 4e-117 | none |

Both are single-copy, unambiguous orthologs in every annotated Pezizales
proteome. APN2 is not found in any of the 7 Tuberaceae proteomes (the
depositors of the Tuber MAT locus also used other anchors); whether it is
lost or unannotated there was not tested.

## Gene order: are they next to MAT?

Counted only where the best MAT hit (E<1e-5) and the flank are on the same
contig (34 genomes):

| flank | on MAT's contig | within 20 kb of MAT | median gap | max gap |
|---|---:|---:|---:|---:|
| SLA2 | 34 | 33 | 6.1 kb | 575 kb (one Tuber assembly) |
| APN2 | 34 | 26 | 16.6 kb | 22 kb |

Morchellaceae (26 genomes): SLA2 4-12 kb and APN2 14-22 kb from MAT, the
KY782629 order (APN2 ... SDH - MAT - MBA1 - end4/SLA2). Tuberaceae: SLA2
14-19 kb from MAT in 4 genomes.

In Pezizaceae, Tarzettaceae, Rhizinaceae and Ascobolaceae the best MAT hit is
weak (E 1e-11 to 1e-14) and never on SLA2's contig: either the annotation
missed the small MAT gene (see the annotation-gap lesson) or the hit is an
HMG-box paralog. Those genomes do not inform gene order.

## Reading

SLA2 is the stronger anchor: universal, single copy, reciprocal best hit to
A. fumigatus everywhere, and 33 of 34 within 20 kb of MAT. APN2 is also an
unambiguous ortholog but sits 14-22 kb out, near the 25 kb clustering gap,
and was not found in Tuberaceae. Even the distant A. fumigatus queries find
both at E <= 1e-85 and 1e-117, so borrowing Pezizomycotina flank proteins
would localise them; the Morchella proteins, from a Pezizales locus deposit,
do so at higher identity (SLA2: 94-99% in Morchellaceae, 58-87% in the other families).
