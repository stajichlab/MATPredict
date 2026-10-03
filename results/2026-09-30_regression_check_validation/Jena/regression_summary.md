# Regression check: curation-umbelopsis a2fe1b4 vs PR #9 8eff88e -- Jena

- Genomes compared: 64; loci rows: 251; loci with a change: 96 (16 touch a call, 80 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 137 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| classifier_shift | 14 |
| idiomorph_changed | 1 |
| span_changed | 2 |

## Changes to calls

- **CBS117697** Mucoromycota:MAT scaffold_4 [span_changed]: called Minus/high -> called Minus/high; margin 57.2 -> 54.1; best score 116.9 -> 112.8
- **CBS156_58** Mucoromycota:MAT scaffold_213 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 103.6 -> 98.4; best score 158.5 -> 154.2
- **CBS206_69** Mucoromycota:MAT scaffold_332 [idiomorph_changed,classifier_shift]: called undetermined/low -> called Minus/low; margin 21.1 -> 30.1; best score 50.8 -> 59.1
- **CBS210_80** Mucoromycota:MAT scaffold_72 [classifier_shift]: called Minus/high -> called Minus/high; margin 86.6 -> 93.3; best score 158.4 -> 162.1
- **CBS223_63** Mucoromycota:MAT scaffold_57 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 84.1 -> 89.4; best score 157.4 -> 159.5
- **CBS251_35** Mucoromycota:MAT scaffold_58 [span_changed]: called Plus/high -> called Plus/high; margin 249.2 -> 248.4; best score 305.9 -> 304.6
- **CBS329_73** Mucoromycota:MAT scaffold_1641 [classifier_shift]: called Minus/high -> called Minus/high; margin 74.0 -> 79.9; best score 151.8 -> 156.5
- **CBS336_62** Mucoromycota:MAT scaffold_4123 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 98.3 -> 108.9; best score 176.4 -> 186.6
- **CBS336_62** Mucoromycota:MAT scaffold_1604 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 91.5 -> 103.0; best score 169.3 -> 180.3
- **CBS338_71** Mucoromycota:MAT scaffold_446 [classifier_shift]: called Minus/high -> called Minus/high; margin 103.7 -> 98.7; best score 148.2 -> 143.8
- **CBS372_39** Mucoromycota:MAT scaffold_532 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 226.7 -> 232.2; best score 294.1 -> 292.8
- **CBS526_68** Mucoromycota:MAT scaffold_680 [classifier_shift]: called Minus/high -> called Minus/high; margin 99.9 -> 106.5; best score 135.5 -> 139.7
- **CBS540_78** Mucoromycota:MAT scaffold_477 [classifier_shift]: called Minus/high -> called Minus/high; margin 95.7 -> 113.7; best score 139.8 -> 149.6
- **CBS541_78** Mucoromycota:MAT scaffold_11 [classifier_shift]: called Plus/high -> called Plus/high; margin 299.1 -> 306.2; best score 357.9 -> 355.1
- **EMLQT1** Mucoromycota:MAT scaffold_19 [classifier_shift]: called Plus/high -> called Plus/high; margin 297.1 -> 302.3; best score 364.8 -> 364.8
- **URM7223** Mucoromycota:MAT scaffold_1367 [classifier_shift]: called Plus/high -> called Plus/high; margin 275.3 -> 280.5; best score 330.2 -> 329.3
