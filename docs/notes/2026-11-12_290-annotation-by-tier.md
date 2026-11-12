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

## What this does and does not show

It supports the v0.9.0 changes at gene level on these 313 genomes: the tier and label rules moved tiers but no gene-level agreement. It says nothing about genomes outside FungiDB, and the disagreements (85 `novel`, 26 `partial` in v0.9.0) have not been looked at one by one; the list is `results/2026-11-12_290_annotation_check/v090_rows.tsv` (HPCC). Plan step 3 (look at the weak-cap trigger loci) and step 5 stay open on #290.
