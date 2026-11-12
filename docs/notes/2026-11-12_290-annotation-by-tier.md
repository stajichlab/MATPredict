# #290: FungiDB annotation check by tier (v0.9.0 run against v0.8.1)

Date 2026-11-12. The annotation is an independent cross-check, not the truth: `novel` means no FungiDB gene overlaps the call (our model may be better, or the annotation may miss the gene), not that the call is wrong.

Inputs: the 313 FungiDB genomes (release 71) run on the frozen v0.8.1 code (`FungiDB_MATpredict`, `--preset v0.8.1`) and on the frozen v0.9.0 code (`7d73be5`, `--preset v0.9.1 --withheld-detail-policy on`, `results/2026-11-09_fdb_v090`); verdicts from `FungiDB_MATpredict/scripts/compare_matpredict_to_gff.py` and `annotation_view.py` run over both; join with `scripts/annotation_by_tier.py`. Results: `results/2026-11-12_290_annotation_check/` (HPCC).

## Agreement by tier (MAT genes in called loci)

| Run | Tier | Genes | Annotated | Partial | Novel |
|---|---|---|---|---|---|
| v0.8.1 | high | 905 | 833 | 15 | 57 |
| v0.8.1 | medium | 329 | 293 | 10 | 26 |
| v0.9.0 | high | 905 | 832 | 15 | 58 |
| v0.9.0 | medium | 344 | 306 | 11 | 27 |

No `low` locus has a gene in these genomes. Every verdict row joins to a report locus.

## v0.9.0 against v0.8.1

- **No verdict changed** on the 1,234 genes present in both runs; none was lost.
- **15 genes only in v0.9.0**, from calls that are new in v0.9.0: 12 in *Bipolaris*, *Curvularia*, *Alternaria* and *Cochliobolus*-type genomes (the `--preset v0.8.1` fix; GAP1, BGL1, MAT1-1-1, MAT1-1-2, MAT1-2-1; all `annotated`) and 3 *Pneumocystis* `matMi` (the corrected record, #302): one `partial`, two `novel`, because the FungiDB models do not cover that ORF.
- **7 genes change tier**: 5 `high` to `medium` (the weak-flank and single-region caps) and 2 `medium` to `high` (*Pneumocystis*, #302). The weak-fragment cap is not in the v0.9.0 preset, so its trigger loci (plan step 3) are not examined here.

## Plan step 3: the weak-fragment cap's trigger loci (not in any preset)

The repo's own rule (`weak_fragment_cap.evaluate`, unchanged) run over the 193 `Ascomycota:MAT` loci of the v0.9.0 FungiDB run with a called idiomorph matches **3** (a one-off script in `results/2026-11-12_290_annotation_check/cap_triggers.py`, not in the repo):

| Genome | Called | Tier | Strongest core gene (the cap's pick) | Annotation at the called idiomorph's core genes |
|---|---|---|---|---|
| *A. terreus* NIH2624 | MAT1-1 | high | MAT1-2-1, polished_single, 44.3 bits: `novel` | MAT1-1-1 `annotated` |
| *C. geophilum* 1.58 | MAT1-1 | medium | MAT1-2-1, polished_single, 45.4 bits: `novel` | MAT1-1-1 `annotated` |
| *S. schenckii* 1099-18 | MAT1-1 | high | MAT1-2-1, polished_single, 37.7 bits: `novel` | MAT1-1-1, MAT1-1-2, MAT1-1-3 all `annotated` |

In all three the cap would pick a thin MAT1-2-1 fragment, a gene of the *other* idiomorph, as the "strongest core gene", while the called idiomorph (MAT1-1) is supported by annotated MAT1-1 genes. So on these genomes the cap would lower the tier of a call the annotation supports. This is three loci in one annotated set, not a measurement of the rule, but it is a reason to read the rule's "strongest core gene" as the strongest gene of the called idiomorph before it is considered for a preset.

## The `novel` and `partial` genes (v0.9.0)

111 genes (85 novel, 26 partial): 107 are `core_MAT` genes (of 689 core genes, 16%) and 4 are flanking genes. The commonest are short ORFs that annotation pipelines often miss: MAT1-2-1 (11), MAT1-1-3 (9), MFalpha (8), MTLA2 (7), RPL39 (7), *Pneumocystis* matMi (6), MAT1-1-1 (5), MAT1-1-2 (4), MAT1-2-4 (4). By family: Ascomycota:MAT 47, Basidiomycota:MAT 23, Ascomycota:MTL 13, Ascomycota:PM 7. They have not been examined genome by genome.

## What this does and does not show

It supports the v0.9.0 changes at gene level on these 313 genomes: the tier and label rules moved tiers but no gene-level agreement. It says nothing about genomes outside FungiDB, and the disagreements have not been looked at one by one; the list is `results/2026-11-12_290_annotation_check/v090_rows.tsv` (HPCC). The one-by-one look at the 111 disagreements and step 5 stay open on #290.
