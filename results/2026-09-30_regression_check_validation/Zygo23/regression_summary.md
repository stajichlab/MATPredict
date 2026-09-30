# Regression check: curation-umbelopsis a2fe1b4 vs PR #9 8eff88e -- Zygo23

- Genomes compared: 46; loci rows: 363; loci with a change: 70 (2 touch a call, 68 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 89 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| span_changed | 2 |

## Changes to calls

- **contig/Chaetocladium_brefeldii_RSA_1136-** Mucoromycota:MAT contig_1326 [span_changed]: called Minus/high -> called Minus/high; margin 33.0 -> 35.5; best score 89.3 -> 91.2
- **contig/Cunninghamella_polymorpha_NRRL_1395** Mucoromycota:MAT contig_634 [span_changed]: called Minus/high -> called Minus/high; margin 39.4 -> 40.0; best score 89.9 -> 89.7
