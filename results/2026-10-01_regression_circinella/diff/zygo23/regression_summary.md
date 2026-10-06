# Regression check: regression run-f9ad092 vs run-a8863f3 -- zygo23

- Genomes compared: 46; loci rows: 382; loci with a change: 57 (24 touch a call, 33 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 81 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| classifier_shift | 16 |
| core_model_changed | 15 |
| span_changed | 2 |

## Changes to calls

- **contig/Absidia_cuneospora_RSA_623_Plus** Mucoromycota:MAT contig_632 [classifier_shift]: called Plus/high -> called Plus/high; margin 174.2 -> 183.1; best score 229.7 -> 242.7
- **contig/Actinomucor_sp._NRRL_A-23671** Mucoromycota:MAT contig_636 [classifier_shift]: called Plus/high -> called Plus/high; margin 335.7 -> 345.4; best score 383.3 -> 395.2
- **contig/Circinomucor_circinelloides_NRRL_22899** Mucoromycota:MAT contig_986 [core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 339.0 -> 346.6; best score 384.8 -> 396.0; models sexM:207bp/28.99%->216bp/30.56%
- **contig/Cunninghamella_bertholletiae_NRRL_1376** Mucoromycota:MAT contig_988 [core_model_changed]: called Minus/high -> called Minus/high; margin 62.9 -> 65.0; best score 119.7 -> 120.4; models sexM:318bp/38.61%->333bp/36.45%
- **contig/Cunninghamella_echinulata_NRRL_1386** Mucoromycota:MAT contig_1220 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 201.9 -> 219.7; best score 272.4 -> 292.6
- **contig/Cunninghamella_echinulata_RSA_2017_Plus** Mucoromycota:MAT contig_527 [classifier_shift]: called Plus/high -> called Plus/high; margin 206.8 -> 222.2; best score 276.5 -> 296.4
- **contig/Cunninghamella_japonica_NRRL_2463** Mucoromycota:MAT contig_1562 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 201.8 -> 220.4; best score 271.0 -> 291.9
- **contig/Mucor_alternans_NRRL_A-15142** Mucoromycota:MAT contig_739 [core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 339.0 -> 346.6; best score 384.8 -> 396.0; models sexM:207bp/28.99%->216bp/30.56%
- **contig/Mucor_circinelloides_var._griseocyanus_NRRL_1416** Mucoromycota:MAT contig_606 [classifier_shift]: called Plus/high -> called Plus/high; margin 322.4 -> 331.4; best score 371.8 -> 386.0
- **contig/Mucor_sp._NRRL_A-14906** Mucoromycota:MAT contig_769 [core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 339.0 -> 346.6; best score 384.8 -> 396.0; models sexM:207bp/28.99%->216bp/30.56%
- **contig/Mucor_sp._NRRL_A-21230** Mucoromycota:MAT contig_747 [core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 339.0 -> 346.6; best score 384.8 -> 396.0; models sexM:207bp/28.99%->216bp/30.56%
- **contig/Mucor_sp._NRRL_A-21232** Mucoromycota:MAT contig_796 [core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 339.0 -> 346.6; best score 384.8 -> 396.0; models sexM:207bp/28.99%->216bp/30.56%
- **contig/Mucor_sp._NRRL_A-21236** Mucoromycota:MAT contig_718 [core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 339.0 -> 347.6; best score 384.8 -> 396.0; models sexM:207bp/28.99%->204bp/30.88%
- **contig/Mucor_sp._NRRL_A-25970** Mucoromycota:MAT contig_695 [core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 339.0 -> 346.6; best score 384.8 -> 396.0; models sexM:207bp/28.99%->216bp/30.56%
- **contig/Phycomyces_blakesleeanus_Burgeff_NRRL_1465** Mucoromycota:MAT contig_509 [classifier_shift]: called Plus/high -> called Plus/high; margin 240.7 -> 230.7; best score 311.6 -> 305.4
- **contig/Phycomyces_blakesleeanus_Burgeff_var._piloboloides_NRRL_2566** Mucoromycota:MAT contig_934 [classifier_shift]: called Plus/high -> called Plus/high; margin 240.7 -> 230.7; best score 311.6 -> 305.4
- **contig/Phycomyces_blakesleeanus_NRRL_1554** Mucoromycota:MAT contig_552 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 240.7 -> 230.7; best score 311.6 -> 305.4
- **scaffold/Circinomucor_circinelloides_NRRL_22899** Mucoromycota:MAT scaffold_511 [core_model_changed]: called Plus/high -> called Plus/high; margin 87.8 -> 89.3; best score 133.6 -> 137.9; models sexM:207bp/28.99%->216bp/30.56%
- **scaffold/Mucor_alternans_NRRL_A-15142** Mucoromycota:MAT scaffold_519 [core_model_changed]: called Plus/high -> called Plus/high; margin 87.8 -> 89.3; best score 133.6 -> 137.9; models sexM:207bp/28.99%->216bp/30.56%
- **scaffold/Mucor_sp._NRRL_A-14906** Mucoromycota:MAT scaffold_525 [core_model_changed]: called Plus/high -> called Plus/high; margin 87.8 -> 89.3; best score 133.6 -> 137.9; models sexM:207bp/28.99%->216bp/30.56%
- **scaffold/Mucor_sp._NRRL_A-21230** Mucoromycota:MAT scaffold_522 [core_model_changed]: called Plus/high -> called Plus/high; margin 87.8 -> 89.3; best score 133.6 -> 137.9; models sexM:207bp/28.99%->216bp/30.56%
- **scaffold/Mucor_sp._NRRL_A-21232** Mucoromycota:MAT scaffold_505 [core_model_changed]: called Plus/high -> called Plus/high; margin 87.8 -> 89.3; best score 133.6 -> 137.9; models sexM:207bp/28.99%->216bp/30.56%
- **scaffold/Mucor_sp._NRRL_A-21236** Mucoromycota:MAT scaffold_485 [core_model_changed]: called Plus/high -> called Plus/high; margin 87.8 -> 84.8; best score 133.6 -> 131.5; models sexM:207bp/28.99%->204bp/30.88%
- **scaffold/Mucor_sp._NRRL_A-25970** Mucoromycota:MAT scaffold_496 [core_model_changed]: called Plus/high -> called Plus/high; margin 87.8 -> 89.3; best score 133.6 -> 137.9; models sexM:207bp/28.99%->216bp/30.56%
