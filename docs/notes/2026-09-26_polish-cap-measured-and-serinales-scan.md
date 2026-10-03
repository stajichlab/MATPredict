# Per-family polish cap: measured; Serinales-wide scan on the flank roster

2026-09-26. Results under `results/2026-09-26_polish_cap/` and
`results/2026-09-26_serinales_all_882aa01/`.

## 1. The cap, measured on a real run

Same code (`f7b9773`), same 567 genomes (five slow panels plus 100
*Saccharomyces*), `--max-polished-clusters-per-family` off vs 6. Rank before
polishing: distinct live genes, then best identity, then hit count
(`compare.py`, `compare_output.txt`).

| panel | calls off | calls cap6 | median s off | median s cap6 | total h off | total h cap6 |
|---|---:|---:|---:|---:|---:|---:|
| dipodascomycetes_etc | 84 | 84 | 1771 | 842 | 44.6 | 21.7 |
| dothideomycetes | 91 | 91 | 620 | 283 | 18.5 | 8.6 |
| orbiliomycetes | 29 | 29 | 790 | 548 | 15.8 | 10.6 |
| pezizomycotina_uncurated | 88 | 88 | 587 | 413 | 15.9 | 10.8 |
| saccharomyces_sample | 210 | 210 | 22 | 25 | 0.6 | 0.7 |
| saccharomycetes_other | 111 | 107 | 1380 | 766 | 39.8 | 22.3 |
| **total (561 genomes with both)** | **613** | **609** | **687** | **419** | **135.3** | **74.6** |

* Compute falls 45% (135.3 -> 74.6 h). The median genome is 39% faster.
* 6 genomes hit the 60-min timeout uncapped; none did with the cap.
* 4 calls are lost, all in the phylum-fallback panel `saccharomycetes_other`:
  * GCA_010994365.1 (*Saccharomycopsis schoenii*): two MTL A `idiomorph_gene_only`
    calls (2 genes, 51%), ranked 8 and 9 of 14 MTL clusters. The genotype
    changes from a+alpha (A high x2, alpha medium) to alpha (medium only).
  * GCA_029290875.1, GCA_030564765.1 (*Saccharomycopsis*): one Ascomycota:MAT
    MAT1-1 medium partial call each -- a Pezizomycotina family called in a
    yeast by phylum fallback. Likely not a true MAT1-1 locus (not checked).

## 2. Rank rule: genes-first stays

`rank_simulation.py` replays a cap on an uncapped run: a call is lost when no
kept cluster of its family overlaps it. On `f7b9773` it reproduces the real
cap6 losses exactly (4 calls, same genomes). It ignores knock-on effects
between clusters, so it is an estimate.

| run (uncapped) | genomes | calls | genes-first N=6 lost | identity-first N=6 lost |
|---|---:|---:|---:|---:|
| cap panels (`f7b9773`), in sample | 561 | 613 | 4 | 0 |
| Serinales-wide (`882aa01`) | 2,368 | 2,647 | 3 | 3 |
| Mucoromycota (`634dda4`) | 293 | 242 | 0 | 2 (Minus, medium) |
| Mortierellomycota / Kickxellomycota | 290 | 13 | 0 | 0 |

Identity-first wins only on the data it was picked from; out of sample it
loses 2 Mucoromycota Minus calls that genes-first keeps. No change to the rank.

The 3 Serinales losses (both rules) are all in *Lodderomyces beijingensis*
GCF_963989305.1: 8 near-identical "A high" MTLa1+a2 blocks on 7 chromosomes;
the cap drops copies 7-9. Genotype unchanged (see report entry D10).

Per N (genes-first; calls lost / share of cluster gene-load still polished):
N=6 4/613, 47%; N=8 2/613, 54%; N=10 1/613, 61% (cap panels). Serinales-wide:
N=6 3/2,647; N=8 1; N=10 0.

**Recommendation (curator decision):** adopt cap 6 as the default. The only
losses measured are in phylum-fallback genomes with many admitted clusters
and in one multi-copy genome. If the fallback losses matter, use N=10 for
phylum-fallback routing only (projected 1 lost call), not a different rank.

## 3. Serinales-wide scan, flank roster (`882aa01`) vs baseline (`8dec9e1`)

`analysis_vs_8dec9e1.txt`. 2,368 genomes; calls 2,647. Main genotype shifts:

* *Candidozyma auris* (952): none 271 -> 12; 258 now alpha.
* *C. albicans* (134): a 93 -> 17, a+alpha 1 -> 77, alpha 0 -> 34, none 39 -> 5.
  Alpha-only assemblies are collapsed heterozygotes (report entry B4).
* *C. tropicalis* (115): a+alpha 0 -> 87. *C. dubliniensis*, *C. viswanathii*,
  *Metschnikowia pulcherrima*, *L. metapsilosis* show the same a -> a+alpha shift.
* Species with no call before: *Lizanozyma spartinae* 110/111 alpha,
  *L. parapsilosis* 97/98 a, *Fermentozyma sake* 10/10, *Candidozyma haemuli*
  10/10, *Diutina catenulata* 17/18.
* MTLalpha2 is found where expected (C. albicans 110 of 115 alpha loci,
  *C. tropicalis* 93/99) and never in *C. auris* or *C. lusitaniae* (D4).

### Flank-carried calls: 131 of 2,647 (4.9%), not 0.6%

The partial scan (910 genomes) showed 7 of 1,096. The full scan shows 131
calls with no modelled core gene (`flank_carried_audit.py`, `.txt`):

| genome also has a core-modelled call? | core hit vs PAP1-OBP1-PIK1 span +-3 kb | calls | main species |
|---|---|---:|---|
| yes | outside | 52 | *D. hansenii* 46, other *Debaryomyces* 5 |
| no | inside | 45 | *Diutina catenulata* 17, *D. rugosa* 4 |
| no | outside | 32 | *L. elongisporus* 5, *Yamadazyma* 10 |
| yes | inside | 2 | |

* ***Debaryomyces*: the flank block is not at the MTL locus.** In CBS767
  (NC_006047.2) PAP1-OBP1-PIK1 is at 0.88 Mb and the a1/a2/alpha1 locus at
  1.59 Mb. The call at the block rests on a 29% MTLA2 fragment 8 kb away. It
  turned 42 `homothallic` genotypes into `a+homothallic`. Report entry D9.
* *L. elongisporus* has no MAT genes (D6); its 5 calls are all "outside".
* The earlier recommended rule handles all of this: core hit outside the flank
  span -> withhold (84 calls, incl. all 52 *Debaryomyces*); inside -> keep at
  `low`, `partial_locus`, `idiomorph_unmodelled: true` (47 calls, mostly
  *Diutina*, which has no curated reference). Not implemented; curator decision.
