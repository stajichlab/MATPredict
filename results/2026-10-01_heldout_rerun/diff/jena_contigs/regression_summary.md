# Regression check: Held-out rerun: 076afe4 vs 52b3ff9 -- jena_contigs

- Genomes compared: 63; loci rows: 243; loci with a change: 115 (27 touch a call, 88 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 113 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| call_gained | 3 |
| call_lost | 6 |
| classifier_shift | 16 |
| gene_set_changed | 2 |
| idiomorph_changed | 1 |
| span_changed | 3 |

## Changes to calls

- **CBS117697** Mucoromycota:MAT contig_6 [span_changed]: called Minus/high -> called Minus/high; margin 56.9 -> 54.1; best score 114.2 -> 112.8
- **CBS169_57** Mucoromycota:MAT contig_4 [call_lost]: called Minus/high -> withheld:paralog_class Minus/; margin 34.1 -> 31.7; best score 81.4 -> 78.3; models no_gene_evidence_in_candidate
- **CBS206_69** Mucoromycota:MAT contig_415 [idiomorph_changed,classifier_shift]: called undetermined/low -> called Minus/low; margin 21.1 -> 30.1; best score 51.1 -> 59.1
- **CBS210_80** Mucoromycota:MAT contig_37 [call_lost]: called Minus/medium -> withheld:paralog_class Minus/; margin 34.1 -> 31.7; best score 81.4 -> 78.3; models no_gene_evidence_in_candidate
- **CBS221_71** Mucoromycota:MAT contig_1 [call_lost]: called Minus/high -> withheld:paralog_class Minus/; margin 35.4 -> 32.9; best score 80.0 -> 76.8; models no_gene_evidence_in_candidate
- **CBS221_71** Mucoromycota:MAT contig_624 [call_gained]: withheld:below_fraction_floor Minus/ -> called Minus/medium; margin 102.2 -> 102.4; best score 168.2 -> 167.2; models no_gene_evidence_in_baseline
- **CBS223_63** Mucoromycota:MAT contig_8 [call_lost]: called Minus/high -> withheld:paralog_class Minus/; margin 33.6 -> 31.3; best score 80.9 -> 77.9; models no_gene_evidence_in_candidate
- **CBS223_63** Mucoromycota:MAT contig_62 [call_gained]: withheld:below_fraction_floor Minus/ -> called Minus/medium; margin 84.6 -> 89.4; best score 156.4 -> 159.5; models no_gene_evidence_in_baseline
- **CBS230_35** Mucoromycota:MAT contig_1186 [classifier_shift]: called Minus/high -> called Minus/high; margin 109.0 -> 103.5; best score 189.1 -> 187.0
- **CBS251_35** Mucoromycota:MAT contig_364 [call_lost,span_changed,gene_set_changed]: called Minus/high -> withheld:paralog_class Minus/; margin 30.5 -> 28.9; best score 76.9 -> 74.4; models no_gene_evidence_in_candidate; genes algA,sexM,tptA -> algA,glrA,sexM,tptA
- **CBS251_35** Mucoromycota:MAT contig_64 [span_changed]: called Plus/high -> called Plus/high; margin 246.5 -> 248.4; best score 303.3 -> 304.6
- **CBS293_63** Mucoromycota:MAT contig_244 [classifier_shift]: called Plus/high -> called Plus/high; margin 309.0 -> 320.6; best score 367.4 -> 372.1
- **CBS329_73** Mucoromycota:MAT contig_1661 [classifier_shift]: called Minus/high -> called Minus/high; margin 73.3 -> 79.9; best score 152.7 -> 156.5
- **CBS336_62** Mucoromycota:MAT contig_4178 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 100.8 -> 108.9; best score 177.7 -> 186.6
- **CBS336_62** Mucoromycota:MAT contig_1643 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 95.0 -> 103.0; best score 170.6 -> 180.3
- **CBS417_77** Mucoromycota:MAT contig_413 [classifier_shift]: called Plus/high -> called Plus/high; margin 327.6 -> 338.8; best score 380.3 -> 384.6
- **CBS421_70** Mucoromycota:MAT contig_3683 [classifier_shift]: called Minus/high -> called Minus/high; margin 126.1 -> 118.3; best score 196.6 -> 188.2
- **CBS526_68** Mucoromycota:MAT contig_739 [classifier_shift]: called Minus/high -> called Minus/high; margin 101.5 -> 106.5; best score 135.7 -> 139.7
- **CBS540_78** Mucoromycota:MAT contig_539 [classifier_shift]: called Minus/high -> called Minus/high; margin 94.9 -> 113.7; best score 139.6 -> 149.6
- **CBS541_78** Mucoromycota:MAT contig_14 [classifier_shift]: called Plus/high -> called Plus/high; margin 295.6 -> 306.2; best score 351.1 -> 355.1
- **CBS576_66** Mucoromycota:MAT contig_408 [classifier_shift]: called Plus/high -> called Plus/high; margin 299.2 -> 304.7; best score 371.6 -> 374.5
- **CBS762_74** Mucoromycota:MAT contig_821 [classifier_shift]: called Minus/high -> called Minus/high; margin 109.4 -> 102.9; best score 146.8 -> 136.4
- **CBS763_74** Mucoromycota:MAT contig_18 [call_lost,gene_set_changed]: called Minus/high -> withheld:paralog_class Minus/; margin 28.5 -> 27.8; best score 84.2 -> 83.3; models no_gene_evidence_in_candidate; genes algA,glrA,rnhA,sexM -> algA,glrA,rnhA,sexM,tptA
- **CBS763_74** Mucoromycota:MAT contig_1048 [call_gained]: withheld:below_fraction_floor Plus/ -> called Plus/medium; margin 255.9 -> 254.5; best score 309.6 -> 309.1; models no_gene_evidence_in_baseline
- **CBS816_70** Mucoromycota:MAT contig_324 [classifier_shift]: called Plus/high -> called Plus/high; margin 177.1 -> 185.9; best score 238.4 -> 243.7
- **EMLQT1** Mucoromycota:MAT contig_35 [classifier_shift]: called Plus/high -> called Plus/high; margin 294.9 -> 302.3; best score 362.0 -> 364.8
- **NRZ2022_0058** Mucoromycota:MAT contig_353 [classifier_shift]: called Plus/high -> called Plus/high; margin 308.0 -> 313.6; best score 358.7 -> 360.6
