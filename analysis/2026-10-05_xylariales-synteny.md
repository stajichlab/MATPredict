# Xylariales SLA2-APN2 neighbourhood: synteny, rearrangement and MAT-region loss

Status: first pass (2026-10-05); RNA-seq coverage and the HMG-box tree are included (section 6 for the tree). Handoff: `docs/HANDOFF-xylariales-2026-10-05.md`.
Data and scripts: `results/2026-10-05_xylariales_nc1011_interval/`.

## Question
*Xylaria flabelliformis* NC1011 (GCA_022453505.1) has the conserved MAT flank genes (APC5, SLA2, COX13, APN2) but no MAT gene
between SLA2 and APN2, and its gene order differs slightly from other Sordariomycetes. Is that an inversion or translocation,
and is the MAT gene missing? Robinson & Natvig 2019 (PMID 30557613) found no canonical MAT in 35 Xylariales.

## Methods in brief
- NC1011 neighbourhood: 14 proteins on JAJLYR010000004.1:225-272 kb (KAI0195597.1-KAI0195610.1; in order ARM, FRE, APC5, CIA30,
  SLA2, COX13, APN2, H604, CPN10, H606, GPR1, H608, H609, H610).
- miniprot (pixi build 0.18-r281; `-N 10 --outs 0.5 --outc 0.3`) of those 14 proteins against 414 BFD genomes: all 257 Xylariales,
  98 Amphisphaeriales (sister order), and up to 10 genera in each of 7 outgroup orders (Sordariales, Hypocreales, Glomerellales,
  Magnaporthales, Ophiostomatales, Diaporthales, Coniochaetales). Hits kept at coverage >= 0.5 and positives >= 0.30. Per genome the
  150-kb window with the most distinct genes is taken, and per gene its best hit there. Orders are normalised so SLA2 reads on the plus strand.
- Two descriptors per genome: the order of SLA2, COX13, APN2 (**SAC** = SLA2, APN2, COX13 with APN2 minus and COX13 plus, the order
  seen in outgroups; **SCA** = SLA2, COX13 minus, APN2 plus, the NC1011 order), and which side of SLA2 APC5 lies on.
- Cross-checks: MATPredict v0.6.0 calls (`results/2026-10-03_ascomycota_v060/loci.tsv`) and the paper's supplementary table
  (`docs/notes/2026-10-05_robinson-natvig-2019-xylariales.tsv`).

## 1. NC1011 interval
Gene models (NCBI): SLA2 KAI0195601.1 (243,984-247,733, partial), COX13 KAI0195602.1 (248,209-249,501, minus), APN2 KAI0195603.1
(249,620-252,216). Gaps: SLA2-COX13 476 bp, COX13-APN2 119 bp, APN2 to the next gene 784 bp, previous gene to SLA2 330 bp.
- None of the 14 annotated proteins in 225-272 kb has an HMG-box, MATalpha_HMGbox or MAT1-1-2 domain (Pfam; E 1e-3).
- None of 139 six-frame ORFs (>= 150 nt) hits those HMMs (E 1e-2).
- tBLASTn of 105 database MAT proteins: best 28 bits (E 0.005), fragments of 20-67 aa. Real Xylariales MAT hits were 39-47 bits.
- RNA-seq (SRR8861595, 70.0 M read pairs, hisat2, 94.99% aligned; unstranded; `work/rna/coverage_summary.tsv`):
  - SLA2 mean 663x, COX13 2,543x (median 708x), APN2 188x, all expressed. The 475-bp SLA2-COX13 gap is transcribed (mean 410x, up to 1,060x), so the SLA2 model,
    partial at its 3' end, is probably truncated and the signal is its 3' end, not a separate gene. The 118-bp COX13-APN2 gap has zero coverage.
  - The 783-bp gap between APN2 and H604 is covered at a flat 57x (20-85x) with no annotated gene and no MAT-domain or tBLASTn hit. It is a stretch
    to inspect (UTR overlap or a small unannotated transcript) and is where a displaced MAT site would lie under the far-side model, but nothing in it
    resembles a MAT gene.
  - Elsewhere in 225-272 kb the only unannotated stretches of 100 bp or more at 10x or higher are next to annotated genes (their UTRs).
So there is no room for, and no sign of, a MAT gene between SLA2 and APN2 in NC1011.

## 2. Gene-order states across 414 genomes
| order | genomes | SAC | SCA | other order | incomplete (< 3 of SLA2, COX13, APN2) | APC5 before / after SLA2 |
|---|---|---|---|---|---|---|
| Xylariales | 257 | 89 | 134 | 12 | 22 | 237 / 11 |
| Amphisphaeriales | 98 | 18 | 14 | 40 | 26 | 95 / 0 |
| Sordariales | 10 | 7 | 0 | 0 | 3 | 0 / 7 |
| Hypocreales | 10 | 5 | 1 | 0 | 4 | 0 / 6 |
| Glomerellales | 9 | 6 | 0 | 1 | 2 | 1 / 6 |
| Diaporthales | 10 | 6 | 1 | 0 | 3 | 0 / 7 |
| Coniochaetales | 1 | 1 | 0 | 0 | 0 | 0 / 1 |
| Magnaporthales | 9 | 0 | 0 | 4 | 5 | 0 / 4 |
| Ophiostomatales | 10 | 0 | 0 | 6 | 4 | not found |
("Before" and "after" are along the normalised strand; APC5 and CIA30 sit on the COX13 side in outgroups and on the opposite side of SLA2
in Xylariales and Amphisphaeriales.)

Xylariales by family (SAC / SCA / other / incomplete): Xylariaceae 15 / 84 / 1 / 7; Hypoxylaceae 54 / 5 / 10 / 8; Diatrypaceae 0 / 44 / 0 / 3;
Microdochiaceae 7 / 0 / 0 / 0; family unassigned 12 / 1 / 1 / 4. Main genera: *Xylaria* 1 SAC / 49 SCA, *Eutypa* 0 / 38, *Nemania* 0 / 9,
*Daldinia* 15 / 2, *Hypoxylon* 16 / 1, *Annulohypoxylon* 16 / 0, *Hypomontagnella* 3 / 2. The state is conserved at genus and family level.
Genome counts overstate independent events: the 257 Xylariales genomes carry 164 species labels, and 41 are *Eutypa lata* strains. By majority state per species,
SAC is 69 species, SCA 73, other 9, incomplete 13.
NC1011 and G536 (same species), *Rosellinia* and *Eutypa lata* are SCA; *M. bolleyi*, *Hypoxylon* CI-4A, *Daldinia* EC12 and *Xylaria* JS573 are SAC.

## 3. Gene-level breakpoints (adjacency of NC1011 neighbours; conserved / genomes with both genes)
| adjacency in NC1011 | outgroups (5 orders) | Amphisphaeriales | Xylariales SAC | Xylariales SCA |
|---|---|---|---|---|
| ARM-FRE | 1/1 | 73/91 (80%) | 84/85 (98%) | 108/108 |
| FRE-APC5 | 0/1 | 1/93 (1%) | 87/88 (98%) | 122/122 |
| APC5-CIA30 | 9/23 (39%) | 95/95 | 88/89 (98%) | 132/133 (99%) |
| **CIA30-SLA2** | **0/20 (0%)** | 95/95 | 89/89 | 133/133 |
| **SLA2-COX13** | **1/28 (3%)** | **0/72** | **0/89** | **86/134 (64%)** |
| COX13-APN2 | 31/31 | 74/74 | 89/89 | 134/134 |
| APN2-H604 | n/a | n/a | 0/2 | 14/14 |
| H604-CPN10 ... H608-H609 | n/a | n/a | n/a | 11/11 to 16/16 |
| H609-H610 | 4/22 (18%) | 30/71 (42%) | 73/86 (84%) | 71/115 (61%) |
- **COX13-APN2 is an invariant unit** in every group, so any event moved or inverted the unit as a whole.
- **Junction 1, CIA30|SLA2 (derived).** Absent in outgroups, present in all Xylariales and Amphisphaeriales: In outgroups APC5 and CIA30 lie on the
  COX13 side of the SLA2-MAT-APN2-COX13 block; in Xylariales and Amphisphaeriales CIA30 is next to SLA2 on the other side. It predates the
  Xylariales-Amphisphaeriales split. Gene order alone cannot tell a translocation from an inversion here.
- **Junction 2, SLA2|COX13 (derived, SCA only).** Absent everywhere except 86 of 134 SCA genomes; in the rest of SCA SLA2 and the unit are not adjacent
  (see the gap classes below). This is the junction that puts the unit, and not the MAT site, next to SLA2.
- Outgroup rows rest on few genomes (distant genomes often lack some of the 14 genes), so treat their percentages as indicative.

## 4. NC1011 against its closest SAC genome (*Xylaria* sp. JS573, the one SAC *Xylaria*)
Read ARM to H610 (JS573 coordinates descend, so its list is reversed):
- NC1011 (SCA): ARM, FRE, APC5, CIA30, SLA2 (244,381-247,733), COX13 (248,567-249,330), APN2 (249,890-252,055), H604, CPN10, H606, GPR1, H608, H609, H610.
- JS573 (SAC): ARM, FRE, APC5, CIA30, SLA2 (1,044,624-1,041,275), 7.5-kb gap, APN2 (1,033,725-1,031,536), COX13 (1,030,919-1,029,937), H609, H610.
RNA-seq weakens the "gain" half of this: three of the five H604-H608 genes show no expression in the SRR8861595 library (H606 0.5x, GPR1 0x, H608 0.8x; H604 243x, CPN10 10,724x), so some of those models may be spurious or condition-specific. This is not a clean inversion of [MAT site, APN2, COX13]: the five genes H604-H608 (~10.5 kb, 253,093-263,646 in NC1011) lie between APN2 and H609 in NC1011 and
have no hit in the JS573 window. The minimal reading is an inversion or translocation of the COX13-APN2 unit relative to SLA2, plus a gain or loss of the
H604-H608 block and loss of the MAT region (JS573's 7.5-kb gap holds the HMG gene the paper placed "between").
Nucleotide breakpoints were not resolved: minimap2 (asm20 and a sensitive preset) on 57-kb windows gave one 704-bp block for NC1011 vs JS573 and none for
*H. argillaceum* (SCA-like, 2.7-kb gap) vs *H. fragiforme* (SAC). Intergenic sequence is too diverged; `breakpoints/summary.md` has the tables.

## 5. Where the MAT-like sequence is
SLA2-APN2 hit-to-hit gap (bp), Xylariales: SAC median 8,898 (21 of 89 under 4 kb); SCA median 2,345 (85 of 134 under 6 kb, i.e. COX13 only; 49 at 6 kb or
more, mostly *Eutypa* 38, *Peroneutypa* 3, with a few *Hypomontagnella*, *Daldinia*, *Diatrype*, *Eutypella*, *Xylaria*, *Alloperoneutypa*).
Correction (2026-10-05, after the dashboard work): the large SCA gaps are mostly a third arrangement, not paralog hits. In 48 SCA genomes another gene lies inside
the SLA2-COX13-APN2 interval: H609 in 44, H610 and H609 in 3, CIA30, H610 and H609 in 1. All 38 *Eutypa lata* genomes with a gap read SLA2, H609, COX13, APN2 with
a gap of about 22.4 kb, so H609 has moved into the interval there; in NC1011 it lies beyond H608.
- **SAC keeps a MAT-sized interval.** 18 of 89 SAC genomes have a v0.6.0 locus within 30 kb of the flank hits, and in 17 of them it lies between SLA2 and the
  APN2/COX13 unit. The paper's two "between" genomes (*M. bolleyi*, *Xylaria* JS573) are both SAC.
- **SCA has no room between SLA2 and the unit.** Of 33 SCA genomes with a v0.6.0 locus at the block, 20 span the unit itself (flank-carried) and 13 extend
  beyond the far side of the unit. None lie on the other side of SLA2. The paper's two matched "adjacent" genomes (*Xylaria* MSU SB201401, *Kretzschmaria
  deusta*) are SCA. That fits a MAT site displaced to the far side of APN2, but the v0.6.0 calls are mostly low or medium confidence and are not independent
  proof of a MAT gene there.
- **A second route to loss exists in SAC.** 21 SAC genomes have a SLA2-APN2 gap under 4 kb (*Daldinia eschscholtzii* group, *D. caldariorum*, *D. bambusicola*,
  some *Hypoxylon*), with no linked HMG gene per the paper for the ones it covers: the MAT region deleted without a reorder.
- The miniprot HMG search in the vicinity (11 references) found only 4 weak hits (positives < 0.5, all SCA) across 223 SAC and SCA genomes. It is too
  insensitive to count as evidence of absence (it also misses the divergent HMG genes the paper reports in *M. bolleyi*); the HMG tree is the better test.
- Paper agreement: "between" 2/2 SAC; "adjacent" 2/2 SCA (1 more incomplete); "unlinked" is mixed (SAC 3, SCA 1, other 3, incomplete 2).

## 6. HMG-box tree placement
Tree: 169 HMG-box domains (HMG_box PF00505 or HMG_box_2 PF09011, plus ~35 aa) from 19 of 31 proteomes (12 had no NCBI proteins), including 91 from Xylariales and
Amphisphaeriales-type genomes of the paper's set and ours; mafft L-INS-i, trimAl -gt 0.3, RAxML-NG LG+G4, 200 bootstraps (`tree/`). Labelled references are only 9 database
MAT1-2-1/MAT1-1-3-type proteins, NCU03481 and fmf-1 (*N. crassa*). Placement: `tree/hmg_placement.tsv` (`08_hmg_placement.py`).
- **KAI0192626.1** (the HMG protein v0.6.0 called MAT1-2 in NC1011): nearest labelled reference is NCU03481 at 1.02 substitutions/site; the nearest MAT reference is 2.22.
  Its smallest clade with a labelled reference has 13 tips and holds NCU03481 only. Bootstrap support for that clade is 14%, so this is a nearest-neighbour result, not a resolved placement.
  It agrees with the expectation that the called gene is an NCU03481-type regulator and not MAT1-2-1.
- Across the 91 Xylariales HMG proteins, 10 are nearest to NCU03481 and 12 to fmf-1. The 69 "nearest to a MAT reference" are uninformative: the median distance to any
  labelled reference is 3.6 substitutions/site (maximum 12.0), i.e. saturated, and most of these are other HMG families (SOX-like and so on).
- **The tree lacks the power to define the MAT1-2-1 clade.** The 9 MAT references are not monophyletic (their smallest common clade has 62 tips, support 21%, and contains NCU03481 and
  fmf-1), unlike the paper's trees, which had many more references and consistently separated MAT1-2-1. A usable test needs more MAT1-2-1/MAT1-1-3 references, NCU03481 and
  fmf-1 orthologs from several Sordariomycetes (by reciprocal best hit), the paper's TreeBase alignment (S23036) as a seed, and the 12 missing proteomes.

## Conclusions
1. The SLA2-APN2 neighbourhood was rearranged twice relative to outgroups (junctions 1 and 2), and the second rearrangement is shared by Xylariaceae and Diatrypaceae
   (NC1011, G536, *Rosellinia*, *Eutypa lata*). So NC1011's order is lineage-wide, not a strain oddity.
2. In that lineage the MAT region is gone from between SLA2 and APN2. Where a MAT-like locus is called it is mostly beyond APN2, which is where a displaced MAT site
   would be.
3. The gene order of the unit COX13-APN2 never breaks, so "inversion or translocation of that unit" is the simplest description; the data cannot separate the two.
4. The earlier note that SAC is "the canonical block inverted" was wrong: relative to outgroups SAC differs by junction 1 only, and SCA adds junction 2.

## Limitations
- Gene-level only; breakpoints are bounded by gene ends, not located to bp. Draft assemblies can mimic a rearrangement; 22 Xylariales genomes were incomplete.
- Window choice, coverage and identity cut-offs are heuristics; paralogs can capture a query. Strain sampling is uneven (41 *Eutypa lata* genomes), so species-level counts are given beside genome counts.
- Outgroup sampling is small and several outgroup orders lack some genes.
- No statistics; counts are descriptive. Calls from v0.6.0 are not validated for these genomes.

## Pending
- A better-powered HMG tree (more references, reciprocal-best-hit NCU03481/fmf-1 orthologs, the 12 proteomes without NCBI annotation, the paper's alignment as seed).
- Closer SAC/SCA pairs for bp-level breakpoints (a protein-level or LASTZ alignment, or long-read assembly of the NC1011/G536 block).
