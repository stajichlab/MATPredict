# Regression check: regression run-9b15d43 vs run-73f56b8 -- zygo23

- Genomes compared: 46; loci rows: 382; loci with a change: 98 (3 touch a call, 95 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 11 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| classifier_shift | 3 |

## Changes to calls

- **scaffold/Cunninghamella_elegans_NRRL_1390** Mucoromycota:MAT scaffold_20 [classifier_shift]: called Minus/high -> called Minus/high; margin 49.6 -> 34.5; best score 106.0 -> 89.5
- **scaffold/Cunninghamella_elegans_NRRL_1393** Mucoromycota:MAT scaffold_446 [classifier_shift]: called Minus/high -> called Minus/high; margin 46.2 -> 34.5; best score 105.0 -> 89.5
- **scaffold/Cunninghamella_elegans_NRRL_2310** Mucoromycota:MAT scaffold_544 [classifier_shift]: called Minus/high -> called Minus/high; margin 46.2 -> 34.5; best score 105.0 -> 89.5
