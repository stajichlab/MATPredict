# Morchella references, SLA2/APN2 in MATtub, and the first homothallic candidates

Curation of 2026-09-25 (commits 81ea22a, ca619b6), measured by re-running the
pilot genomes from frozen worktree `run-d689647`. Comparisons:
`results/2026-09-25_{pezizomycetes,pichiomycetes}_d689647/compare_vs_2026-09-24.txt`
(script `results/compare_panels.py`).

## What changed in the database

* `1174673_ypl6-1_MATtub_MAT1-2` and `1174673_ypl6-3_MATtub_MAT1-1`, Morchella
  importuna, from the Chai et al. 2017 locus deposits KY782629.1 / KY782630.1
  (doi:10.1007/s11557-017-1309-x; no PMID, curator ruling recorded). 7/7 genes
  validate at 100% identity and coverage.
* SLA2 and APN2 added to the MATtub roster as optional flanks.
* `4959_cbs767_MTL_A` gains MTLalpha1 (XP_460134.1), record_version 2,
  `system: homothallic` as a candidate from gene content.

## A runner defect found on the way, fixed in d689647

The first re-run showed no change at all. `detect` reads `$MATPREDICT_DB_ROOT`,
else `<working dir>/db`, and a SLURM job runs in its submission directory -- the
main checkout, on an older branch. The worktree's code had searched the main
checkout's database, without the new records. Both runners now pin the
database to the tree the code comes from. Earlier pilots and the holdout were
unaffected in practice: the only database difference between the two trees
was locus-tag metadata on two records, which changes no reference protein.

## Pezizomycetes, 172 genomes

| family | before | after | new calls adjacent to SLA2/APN2 |
|---|---:|---:|---:|
| Discinaceae | 1/72 | **41/72** | 39 |
| Morchellaceae | 14/60 | **43/60** | 43 |
| Tuberaceae | 12/12 | 12/12 | 5 |
| Pyronemataceae | 0/9 | 2/9 | 2 |
| Pezizaceae, Rhizinaceae, Tarzettaceae, Ascodesmidaceae, Ascobolaceae | 0/19 | 0/19 | - |
| **total** | **27** | **98** | 89 of 98 |

* No call was lost: all 27 earlier calls are still overlapped by a new call on
  the same contig with the same idiomorph. 19 grew outward to take in a flank.
* 79 of 98 calls now include SLA2. Adjacency is to each genome's SLA2/APN2
  ortholog located independently of the pipeline
  (`results/2026-09-24_bar_synteny/`). Of the 9 non-adjacent calls, 7 are the
  existing Tuber calls, where SLA2 sits on another contig.
* Most new calls are `medium` (high-confidence 27 -> 29). Idiomorphs: MAT1-1 60,
  MAT1-2 38.
* Cost: median wall time 22 s -> 60 s per genome (p90 92 s, max 204 s),
  from polishing the added flanks.

## Pichiomycetes, 95 of 99 genomes

Debaryomycetaceae 2 -> 6 genomes called; the first `alpha` calls appear.
Metschnikowiaceae stays 0/30 and Pichiaceae 12/20, as expected: the only new
MTL reference is D. hansenii.

## Homothallic candidates the current rule does not label

`classify_locus` (idiomorph.py:472) emits `homothallic_candidate` only when BOTH
opposite-idiomorph core genes come from an annotated proteome. The guard exists
because sexM/sexP are both HMG-box genes, and on genome-only runs two partial
alignments of one HMG region looked like two genes (54 of 75 false candidates
in the 44-genus sweep). It therefore cannot fire on any genome-only run.

These re-runs contain loci that look like the real thing, where the two genes
are UNRELATED proteins and so cannot be one region seen twice:

| genome | species | genes at one locus | labelled |
|---|---|---|---|
| GCA_040803765.1 | Hydnotrya cerebriformis | MAT1-2-1 (2,578-3,156) + MAT1-1-1 (5,139-6,075), both exonerate models, ~49% | MAT1-2 |
| GCA_040803455.1 | Hydnotrya variiformis | MAT1-1-1 (8,931-9,867) + MAT1-2-1 (11,840-12,418), both exonerate models, ~49% | MAT1-2 |
| GCA_030564605.1 | Debaryomyces coudertii NRRL Y-7425 | MTLA1 + MTLA2 + MTLalpha1 | alpha |
| GCA_056149625.1 | Debaryomyces hansenii Wch | MTLA1 + MTLA2 + MTLalpha1 | A |
| GCA_030574455.1 | Priceomyces fermenticarens NRRL Y-17321 | MTLA1 + MTLA2 + MTLalpha1 | alpha |
| GCA_030574515.1 | Priceomyces melissophilus NRRL Y-7585 | MTLA2 + MTLalpha1 | alpha |
| GCA_046867635.1 | Schwanniomyces etchellsii CBS 6823 | MTLA2 + MTLalpha1 | alpha |
| GCA_054122755.1 | Debaryomyces hansenii (MAG AOP7_A2_surf_bin_9) | MTLA1 + MTLalpha1 (two loci) | alpha |

MAT1-1-1 (alpha box) and MAT1-2-1 (HMG box) are not homologous, nor are MTLa1
(homeodomain) and MTLalpha1 (alpha box). Dirks et al. 2025 report one
colocalised MAT1-1/MAT1-2 case in Discinaceae needing confirmation; Krassowski
et al. 2019 describe the contiguous a+alpha arrangement in Debaryomyces.
Whether these genomes self-mate is not shown by any of this.

Proposed (curator decision, not implemented): accept polished models in the
homothallic rule when the two genes' `gene_class` or domain differ, keep the
proteome-only guard for same-domain pairs such as sexM/sexP, and require the
two models not to overlap.
