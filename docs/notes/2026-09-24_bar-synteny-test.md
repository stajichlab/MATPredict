# Is the modelled-gene bar withholding real loci? A synteny test

Data: the 841 pilot genomes (`results/2026-09-24_pilots/`), every withheld locus
from `suppressed_loci`. Test: `results/2026-09-24_bar_synteny/` (`find_flanks.sh`,
`analyze.py`, `adjacency_*_core.txt`).

## Why a separate test was needed

"Nearly every uncalled genome has a bar-withheld locus" (the pilot note) mixes
two things. Withheld loci include genome-wide background: in Serinales every
cluster was a single homeodomain HSP, never admitted to polishing. And in the
Pezizomycotina `MAT` family, background in CALLED genomes is itself 100%
core+flank and 100% admitted, at ~36% identity, so structure and admission do
not separate a real locus from a paralog cluster there.

Position does. Pezizomycotina and many yeasts keep MAT next to SLA2 and/or
APN2. For each genome the SLA2 and APN2 orthologs were located independently
of the pipeline: the top tblastn hit (bitscore) of the A. fumigatus and
Diaporthales proteins. A locus is ADJACENT when it carries a core MAT gene and
lies within 20 kb of either ortholog on the same contig. The core-gene
requirement matters because SLA2/APN2 are roster genes in some families.

## The test is valid where the linkage holds

| pilot | called genomes with an adjacent call | background loci adjacent |
|---|---:|---:|
| uncurated Pezizomycotina orders | 82/82 (100%) | 2/616 |
| Dipodascomycetes etc. | 47/48 (98%) | 115/5178 (2%) |
| Dothideomycetes | 86/90 (96%) | 7/920 (1%) |
| Orbiliomycetes | 26/27 (96%) | 18/2784 (1%) |
| Saccharomycetes, not Saccharomycetaceae | 51/56 (91%) | 202/6483 (3%) |
| Pezizomycetes | 20/27 (74%) | 0/2 |
| Pichiomycetes | 11/18 (61%) | 60/1889 (3%) |
| Taphrinomycotina | 13/80 (16%) | 60/2257 (3%) |

Called loci are adjacent 91-100% of the time where SLA2/APN2 linkage is
conserved; background loci 0-3%. The test does NOT apply to
Taphrinomycotina (16%) and is weak in Pichiomycetes (61%; the CTG-clade MTL
is flanked by PAP/OBP/PIK) and Pezizomycetes (74%). Results were stable at
10 kb and 50 kb windows except Orbiliales and Taphrinomycotina.

## Result: the bar is the problem in some lineages, the references in others

Genomes with no call, and whether ANY withheld locus is adjacent:

| order | withheld-only genomes | with an adjacent core locus | reading |
|---|---:|---:|---|
| Pezizales (outside Tuberaceae) | 81 | **69 (85%)** | bar withholds the real locus |
| Dipodascales | 46 | **35 (76%)** | bar withholds the real locus |
| Saccharomycodales | 18 | **15 (83%)** | bar withholds a core MAT gene next to SLA2 |
| Pichiales | 8 | 7 | bar |
| Trigonopsidales | 4 | 4 | bar |
| Ascoideales | 15 | 8 | mixed |
| Phaffomycetales | 9 | 5 | mixed |
| Taphrinales | 17 | 6 | test weak here |
| Orbiliales | 40 | 5 (12%) | locus NOT localised: reference gap |
| Serinales | 73 | 3 (4%) | locus NOT localised: reference gap |
| Dothideomycetes | 10 | 2 | mostly not localised |
| uncurated Pezizomycotina | 18 | 3 | mostly not localised |

Saccharomycodales is Hanseniaspora, where MAT gene loss has been proposed.
15 of 18 no-call genomes carry a core MAT hit next to SLA2. That is evidence
against a blanket loss, not a verification of intact genes.

## What lowering the bar would do, versus a synteny rule

A lower bar reports exactly the withheld loci, so its effect can be read off:

| pilot | bar=1: genomes gained | loci | adjacent loci |
|---|---:|---:|---:|
| Pezizomycetes | 41 | 41 | 37 (90%) |
| Dipodascomycetes etc. | 44 | 131 | 77 (59%) |
| Saccharomycetes other | 38 | 125 | 54 (43%) |
| Orbiliomycetes | 10 | 13 | 5 |
| Pichiomycetes | 9 | 34 | 1 |

`bar=0` gains more genomes but reports thousands of background loci (e.g.
5,528 in 51 Dipodascomycetes genomes, 3% adjacent). A global bar change is
therefore not clean outside Pezizales: bar=1 roughly triples the loci per
gained genome in the yeasts. **Synteny is the discriminator**: a withheld locus
with a core gene next to the genome's SLA2/APN2 ortholog is real at a
background rate of 0-3% per candidate locus in the lineages where the test is
valid.

## Options (curator decision)

1. **Synteny rescue tier.** Report a withheld locus, capped (low, its own
   class), when it carries a core gene and lies near the genome's top SLA2 or
   APN2 hit. SLA2/APN2 are already roster genes in `MAT` and `MATsc`; `MATtub`
   would need them added (Morchella flanks include both). Would recover about
   69 Pezizales, 35 Dipodascales, 15 Saccharomycodales genomes in these pilots.
2. **Curation for the lineages where nothing is localised**: Orbiliales,
   Serinales, the residual Dothideomycetes/uncurated Pezizomycotina. No bar
   change helps these.
3. Leave the bar, report bar-only genomes as "withheld candidate present", so a
   sweep does not read them as absence.
