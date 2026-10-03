# Regression check: regression 94d9d37 vs a8863f3 (baseline) -- mucoromycota

- Genomes compared: 80; loci rows: 944; loci with a change: 279 (28 touch a call, 251 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 396 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| call_gained | 3 |
| classifier_shift | 24 |
| core_model_changed | 8 |
| gene_set_changed | 3 |
| span_changed | 10 |

## Changes to calls

- **GCA_000325505.1_RHIrdgD1.0** Mucoromycota:MAT ANKS01000730.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 143.5 -> 128.5; best score 216.4 -> 199.8
- **GCA_000401635.1_Muco_sp_1006Ph_V1** Mucoromycota:MAT KE123975.1 [core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 339.0 -> 346.6; best score 384.8 -> 396.0; models sexM:207bp/28.99%->216bp/30.56%
- **GCA_000696955.1_SynRacB6101-1.0** Mucoromycota:MAT JNDN01000979.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 189.7 -> 182.8; best score 227.7 -> 219.1
- **GCA_000697035.1_RhiStoB9770-1.0** Mucoromycota:MAT JNDS01005206.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 280.2 -> 292.3; best score 329.0 -> 341.1
- **GCA_000697055.1_SakVasB4078-1.0** Mucoromycota:MAT JNDT01002012.1 [span_changed,classifier_shift]: called Minus/high -> called Minus/high; margin 61.5 -> 70.0; best score 141.2 -> 148.2
- **GCA_000697215.1_CunBer175-1.0** Mucoromycota:MAT JNEG01000802.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 207.3 -> 219.6; best score 266.7 -> 284.5
- **GCA_000697315.1_CunBerB7461-1.0** Mucoromycota:MAT JNEL01000814.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 207.3 -> 219.6; best score 266.7 -> 284.5
- **GCA_002105135.1_Synrac1** Mucoromycota:MAT MCGN01000004.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 189.7 -> 182.8; best score 227.7 -> 219.1
- **GCA_013461545.1_ASM1346154v1** Mucoromycota:MAT VAFG01000366.1 [core_model_changed]: called Minus/medium -> called Minus/medium; margin 35.1 -> 38.8; best score 105.7 -> 107.4; models sexM:246bp/41.46%->246bp/41.46%
- **GCA_016758895.1_ASM1675889v1** Mucoromycota:MAT JAEPRA010000016.1 [span_changed,classifier_shift]: called Minus/high -> called Minus/high; margin 100.6 -> 82.5; best score 127.3 -> 110.9
- **GCA_016758965.1_ASM1675896v1** Mucoromycota:MAT JAEPRB010000020.1 [call_gained,span_changed,classifier_shift]: withheld:mat_gene_gate Plus/ -> called Plus/medium; margin 54.2 -> 134.7; best score 86.4 -> 165.8; models no_gene_evidence_in_baseline
- **GCA_021235465.1_Umbelo1** Mucoromycota:MAT JAJSME010000005.1 [span_changed,gene_set_changed]: called Plus/high -> called Plus/high; margin 135.3 -> 138.2; best score 165.6 -> 173.9; genes algA,glrA,sexP,tptA -> algA,glrA,rnhA,sexP,tptA
- **GCA_023630305.1_ASM2363030v1** Mucoromycota:MAT JAMAMS010000120.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 164.2 -> 170.5; best score 195.8 -> 200.6
- **GCA_024139275.1_ASM2413927v1** Mucoromycota:MAT VCJA01002019.1 [core_model_changed,classifier_shift]: called Plus/medium -> called Plus/medium; margin 252.7 -> 244.9; best score 318.4 -> 312.7; models sexM:243bp/30.86%->249bp/30.12%
- **GCA_024139395.1_ASM2413939v1** Mucoromycota:MAT VCJI01001366.1 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 94.0 -> 84.0; best score 166.7 -> 157.5
- **GCA_024139405.1_ASM2413940v1** Mucoromycota:MAT VCJL01001223.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 89.3 -> 80.2; best score 163.8 -> 156.3
- **GCA_025093555.1_Cirumb1** Mucoromycota:MAT JAIWNF010000086.1 [call_gained,span_changed,classifier_shift]: withheld:mat_gene_gate Plus/ -> called Plus/medium; margin 54.3 -> 144.2; best score 88.3 -> 176.1; models no_gene_evidence_in_baseline
- **GCA_025266875.1_Fenlin1** Mucoromycota:MAT JALLLT010000045.1 [span_changed,core_model_changed,classifier_shift]: called Minus/medium -> called Minus/medium; margin 56.4 -> 85.0; best score 121.8 -> 149.3; models sexM:360bp/52.99%->558bp/59.14%
- **GCA_025531835.1_PnitS608_1** Mucoromycota:MAT JAJDOM010000040.1 [core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 241.1 -> 227.2; best score 308.3 -> 297.3; models sexM:219bp/36.99%->237bp/34.18%
- **GCA_025677815.1_Phaart1** Mucoromycota:MAT JAIXMP010000004.1 [call_gained,span_changed,classifier_shift]: withheld:mat_gene_gate undetermined/ -> called Plus/medium; margin 22.3 -> 92.7; best score 48.5 -> 116.9; models no_gene_evidence_in_baseline
- **GCA_037042155.1_ASM3704215v1** Mucoromycota:MAT JAXQJX010000258.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 164.3 -> 170.0; best score 195.6 -> 199.9
- **GCA_054906655.1_Umbelopsis_nana_v._1.0** Mucoromycota:MAT BAAHSD010000001.1 [span_changed,gene_set_changed]: called Minus/high -> called Minus/high; margin 31.2 -> 28.7; best score 50.8 -> 48.4; genes algA,glrA,sexM,tptA -> algA,glrA,rnhA,sexM,tptA
- **GCA_059714295.1_ASM5971429v1** Mucoromycota:MAT JCAOFY010000002.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 154.3 -> 161.0; best score 190.7 -> 197.6
- **GCA_977110945.1_gzUmbRama1** Mucoromycota:MAT CDSBDH010000020.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 137.9 -> 144.4; best score 169.7 -> 182.3
- **GCA_977110975.1_gzUmbVina2** Mucoromycota:MAT CDSBDG010000016.1 [gene_set_changed]: called Minus/high -> called Minus/high; margin 36.8 -> 33.2; best score 58.3 -> 56.0; genes algA,glrA,sexM,sexP,tptA -> algA,glrA,rnhA,sexM,sexP,tptA
- **GCF_025331425.1_Radspe1** Mucoromycota:MAT NW_026251940.1 [core_model_changed,classifier_shift]: called Plus/medium -> called Plus/medium; margin 172.9 -> 166.8; best score 307.7 -> 313.6; models sexM:249bp/55.42%->552bp/41.99%
- **GCF_025528875.1_Mycafr1** Mucoromycota:MAT NW_026515730.1 [core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 239.0 -> 227.1; best score 294.7 -> 284.2; models sexM:237bp/30.77%->216bp/31.94%
- **GCF_025766255.1_Zycmex1** Mucoromycota:MAT NW_026516701.1 [span_changed,core_model_changed,classifier_shift]: called Minus/medium -> called Minus/medium; margin 61.7 -> 133.6; best score 132.6 -> 203.7; models sexM:360bp/58.97%->573bp/100.0%
