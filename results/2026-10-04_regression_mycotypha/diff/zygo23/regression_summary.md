# Regression check: regression 356b9f3 (Mycotypha record + rebuild) vs v0.6.1 (b04a77b) -- zygo23

- Genomes compared: 46; loci rows: 381; loci with a change: 12 (6 touch a call, 6 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 108 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| classifier_shift | 4 |
| span_changed | 2 |

## Changes to calls

- **contig/Absidia_cuneospora_RSA_623_Plus** Mucoromycota:MAT contig_632 [classifier_shift]: called Plus/high -> called Plus/high; margin 183.1 -> 176.8; best score 242.7 -> 238.5
- **contig/Cunninghamella_bertholletiae_NRRL_1376** Mucoromycota:MAT contig_988 [span_changed]: called Minus/high -> called Minus/high; margin 65.0 -> 61.1; best score 120.4 -> 116.7
- **contig/Cunninghamella_bertholletiae_NRRL_1380** Mucoromycota:MAT contig_478 [span_changed]: called Minus/high -> called Minus/high; margin 67.5 -> 63.7; best score 123.5 -> 120.5
- **contig/Cunninghamella_echinulata_NRRL_1386** Mucoromycota:MAT contig_1220 [classifier_shift]: called Plus/high -> called Plus/high; margin 219.7 -> 205.3; best score 292.6 -> 278.9
- **contig/Cunninghamella_echinulata_RSA_2017_Plus** Mucoromycota:MAT contig_527 [classifier_shift]: called Plus/high -> called Plus/high; margin 222.2 -> 208.3; best score 296.4 -> 282.6
- **contig/Cunninghamella_japonica_NRRL_2463** Mucoromycota:MAT contig_1562 [classifier_shift]: called Plus/high -> called Plus/high; margin 220.4 -> 205.9; best score 291.9 -> 278.2
