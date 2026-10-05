# Regression check: regression 84b6b55 (frameshift-aware classification) vs v0.6.1 (b04a77b) -- mucoromycota

- Genomes compared: 80; loci rows: 1007; loci with a change: 190 (12 touch a call, 178 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 496 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| call_gained | 2 |
| classifier_input_changed | 2 |
| classifier_shift | 11 |
| core_model_changed | 1 |
| span_changed | 5 |

## Changes to calls

- **GCA_000697015.1_CunEleB9769-1.0** Mucoromycota:MAT JNDR01001106.1 [span_changed]: called Minus/high -> called Minus/high; margin 67.5 -> 63.7; best score 123.5 -> 120.5
- **GCA_000697035.1_RhiStoB9770-1.0** Mucoromycota:MAT JNDS01005206.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 292.3 -> 285.6; best score 341.1 -> 338.6
- **GCA_000697215.1_CunBer175-1.0** Mucoromycota:MAT JNEG01000802.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 219.6 -> 203.1; best score 284.5 -> 270.1
- **GCA_000697315.1_CunBerB7461-1.0** Mucoromycota:MAT JNEL01000814.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 219.6 -> 203.1; best score 284.5 -> 270.1
- **GCA_001683725.1_ASM168372v1** Mucoromycota:MAT LUGH01000105.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 271.0 -> 223.8; best score 330.3 -> 286.4
- **GCA_024139405.1_ASM2413940v1** Mucoromycota:MAT VCJL01001223.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 80.2 -> 87.9; best score 156.3 -> 162.8
- **GCA_025331435.1_Chocucu1** Mucoromycota:MAT JAIWNI010000021.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 268.2 -> 220.5; best score 328.7 -> 284.4
- **GCA_900175165.2_FCH_5_7** Mucoromycota:MAT FWWN02000384.1 [call_gained,span_changed,classifier_input_changed,classifier_shift]: withheld:modelled_gene_bar Plus/ -> called Plus/medium; margin 72.8 -> 85.5; best score 200.1 -> 195.3; models no_gene_evidence_in_baseline
- **GCF_025331425.1_Radspe1** Mucoromycota:MAT NW_026251940.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 166.8 -> 173.9; best score 313.6 -> 317.6
- **GCF_025528875.1_Mycafr1** Mucoromycota:MAT NW_026515730.1 [core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 227.1 -> 175.5; best score 284.2 -> 233.4; models sexP:525bp/45.71%->741bp/100.0%
- **GCF_025766255.1_Zycmex1** Mucoromycota:MAT NW_026516701.1 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 133.6 -> 127.7; best score 203.7 -> 198.0
- **GCF_041956525.1_Rhipu1** Mucoromycota:MAT NW_027192143.1 [call_gained,span_changed,classifier_input_changed,classifier_shift]: withheld:modelled_gene_bar Plus/ -> called Plus/medium; margin 69.7 -> 79.7; best score 202.5 -> 195.3; models no_gene_evidence_in_baseline
