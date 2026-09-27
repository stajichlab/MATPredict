# Mucoromycota-group scan and the chytrid negative control

2026-09-26, code `634dda4` (frozen worktree `run-634dda4`), results
`results/2026-09-26_early_diverging/` (`summary.txt`, `withheld_profile.txt`).

## Chytrid negative control: no false calls; exhaustive routing is the cost

25 chytrids across 10 orders, default routing. No chytrid MAT locus is known, and
no chytrid family is curated, so every genome falls to `exhaustive` routing
(every family of every phylum).

* **0 calls in 19 finished genomes.** No false positive at the report level.
  Every genome had withheld loci (background hits that fail the bar).
* **6 of 25 hit the 60-minute per-genome timeout**; the 19 that finished took a
  median 1,866 s (max 3,121 s). For a public sweep, any genome outside a
  curated phylum costs 30-60 min and yields nothing.
* Implication for the sweep: route genomes of phyla with no curated family to
  "not searched" (a report that says so) rather than `exhaustive`, or cap
  their runtime. Curator decision; not implemented.

## Mucoromycota: 227/293 genomes called

| family | called | notes |
|---|---:|---|
| Mucoraceae | 69/70 | Plus 36 / Minus 39 loci |
| Rhizopodaceae | 93/111 | Plus-biased (see below) |
| Cunninghamellaceae | 18/22 | |
| Umbelopsidaceae | 13/14 | |
| Backusellaceae | 13/14 | |
| Choanephoraceae | 7/7 | |
| Phycomycetaceae | 5/6 | |
| **Lichtheimiaceae** | **1/30** | 29 withheld-only; 23 of them core+flank, 0 modelled, ~42% identity |
| Syncephalastraceae | 4/9 | 5 withheld-only, core+flank |
| Endogonales | 0/5 | |

* **Lichtheimiaceae is the clearest curation target** (Lichtheimia 17,
  Rhizomucor 5, Zychaea, Circinella, Thermomucor): the locus is found with its
  flanks but nothing models at ~42% identity -- the Pezizales pattern. No
  Lichtheimiaceae record exists.
* **Rhizopus "Plus bias" is largely a LABELLING ERROR, not biology or sampling.**
  (Revised 2026-09-26; the earlier text here called the skew real and the
  curator's first view was "likely sampling bias". Both are withdrawn.)
  33 Mucoromycota calls carry sexM (about 96% identity, 99% coverage) and no
  sexP, yet are labelled Plus. All 33 are Rhizopus: R. arrhizus 19,
  R. delemar 9, R. stolonifer 3, R. microsporus 2. 19 of the 43 R. arrhizus
  Plus calls are of this kind.
  **Where the label comes from:** `idiomorph_candidates`
  (`src/MATPredict/detect/idiomorph.py`) scores each idiomorph by the best
  bitscore of any gene with `present_in_idiomorphs` that is not marked
  `idiomorph_informative: false`. In `db/Mucoromycota/order.yml` the flank
  gene btbA is `present_in_idiomorphs: ["Plus"]`, because the only curated
  record that carries btbA (CBS 346-36) is a Plus strain. btbA is present in
  all 33 of these calls. It is a long protein, so at ~98% identity its
  bitscore (e.g. 871) beats the sexM bitscore (e.g. 378). A flank gene outvotes
  the core MAT gene. Example: GCA_000696915.1, Plus 871 vs Minus 378, while
  the per-contig sexM/sexP resolution names sexM the winner (95.7% vs 32.1%).
  btbA is in 74 Plus calls and 1 Minus call overall.
  The Plus/Minus ratio in Rhizopus cannot be read until this is fixed and the
  scan re-scored. Open: 12 R. arrhizus genomes stay uncalled at ~40% identity
  -- assembly fragmentation or divergent references, not checked.

## Mortierellomycota (6/100) and Kickxellomycota (7/190): leads, not results

Searched with the Mucoromycota family (`--phylum Mucoromycota`); no MAT locus is
described for either phylum, so there is no truth to score against.

* Most genomes have a withheld locus (93/100, 179/190), most often sexM/sexP
  hits co-located with tptA/rnhA/glrA flanks at 37-46% identity.
* The 13 calls carry a modelled sexM with flanks (e.g. Mortierella alpina:
  sexM 43%, glrA 36%, rnhA 26%; Linnemannia gamsii: glrA 66%, rnhA 73%, sexP
  29%). Several Coemansia calls have identical values, consistent with a
  shared locus in closely related strains or with duplicated assemblies.
* The round-2 literature check found no tptA-HMG-rnhA synteny in these phyla.
  These co-locations are therefore either a discovery or HMG paralogs next to
  flank paralogs. **Before believing them:** confirm the flanks are the
  genomes' true tptA/rnhA orthologs (reciprocal best hits) and that the HMG
  gene sits between them, as the Pezizomycotina SLA2/APN2 test did.

## Curation targets from this scan, in order

1. Lichtheimiaceae (30 genomes, 1 called) -- find a Lichtheimia/Rhizomucor
   sexM/sexP locus deposit.
2. Syncephalastraceae (9, 4 called) -- a Syncephalastrum reference.
3. Mortierellomycota / Kickxellomycota -- first a synteny check (above); only
   then decide whether they are a curation target or a discovery project.
4. Endogonales (5, 0 called).
