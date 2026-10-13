# Regression check: regression run-9b15d43 vs run-73f56b8 -- mucoromycota

- Genomes compared: 80; loci rows: 1064; loci with a change: 999 (17 touch a call, 982 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 1 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| confidence_changed | 3 |
| core_model_changed | 2 |
| gene_set_changed | 10 |
| locus_class_changed | 1 |
| span_changed | 17 |

## Changes to calls

- **GCA_000325505.1_RHIrdgD1.0** Mucoromycota:MAT ANKS01000730.1 [span_changed]: called Minus/high -> called Minus/high; margin 124.3 -> 124.3; best score 196.1 -> 196.1
- **GCA_000587855.1_B50** Mucoromycota:MAT KK076501.1 [span_changed]: called Minus/high -> called Minus/high; margin 39.6 -> 39.6; best score 89.8 -> 89.8
- **GCA_000696895.1_MucVelB5328-1.0** Mucoromycota:MAT JNDK01002395.1 [span_changed]: called Minus/high -> called Minus/high; margin 106.8 -> 106.8; best score 146.3 -> 146.3
- **GCA_000696955.1_SynRacB6101-1.0** Mucoromycota:MAT JNDN01000979.1 [confidence_changed,locus_class_changed,span_changed,gene_set_changed]: called Plus/high -> called Plus/medium; margin 182.2 -> 182.2; best score 218.9 -> 218.9; genes algA,rnhA,sexP -> rnhA,sexP
- **GCA_002105135.1_Synrac1** Mucoromycota:MAT MCGN01000004.1 [span_changed]: called Plus/medium -> called Plus/medium; margin 182.2 -> 182.2; best score 218.9 -> 218.9
- **GCA_013461545.1_ASM1346154v1** Mucoromycota:MAT VAFG01000366.1 [span_changed]: called Minus/medium -> called Minus/medium; margin 39.4 -> 39.4; best score 108.0 -> 108.0
- **GCA_016758895.1_ASM1675889v1** Mucoromycota:MAT JAEPRA010000016.1 [span_changed,gene_set_changed,core_model_changed]: called Minus/high -> called Minus/high; margin 80.5 -> 80.5; best score 109.5 -> 109.5; models sexP:removed; genes algA,glrA,rnhA,sexM,tptA -> algA,glrA,sexM,tptA
- **GCA_016758905.1_ASM1675890v1** Mucoromycota:MAT JAEPQZ010000017.1 [confidence_changed,span_changed,gene_set_changed]: called Plus/medium -> called Plus/high; margin 76.0 -> 76.0; best score 113.3 -> 113.3; genes algA,sexM,sexP,tptA -> algA,sexP,tptA
- **GCA_021235465.1_Umbelo1** Mucoromycota:MAT JAJSME010000005.1 [span_changed,gene_set_changed]: called Plus/high -> called Plus/high; margin 138.6 -> 138.6; best score 177.7 -> 177.7; genes algA,glrA,rnhA,sexP,tptA -> algA,glrA,sexP,tptA
- **GCA_024139405.1_ASM2413940v1** Mucoromycota:MAT VCJL01001223.1 [confidence_changed,span_changed,gene_set_changed]: called Minus/high -> called Minus/medium; margin 87.9 -> 87.9; best score 162.8 -> 162.8; genes algA,glrA,sexM,tptA -> algA,glrA,sexM
- **GCA_025266875.1_Fenlin1** Mucoromycota:MAT JALLLT010000045.1 [span_changed]: called Minus/medium -> called Minus/medium; margin 85.9 -> 85.9; best score 149.2 -> 149.2
- **GCA_025528865.1_Pilumb1** Mucoromycota:MAT JAIXMH010000022.1 [span_changed]: called Plus/high -> called Plus/high; margin 148.6 -> 148.6; best score 209.0 -> 209.0
- **GCA_037039815.1_ASM3703981v1** Mucoromycota:MAT JAXQGL010000038.1 [span_changed,gene_set_changed]: called Plus/high -> called Plus/high; margin 76.5 -> 76.5; best score 107.8 -> 107.8; genes algA,glrA,sexM,sexP,tptA -> algA,sexM,sexP,tptA
- **GCA_054906655.1_Umbelopsis_nana_v._1.0** Mucoromycota:MAT BAAHSD010000001.1 [span_changed,gene_set_changed,core_model_changed]: called Minus/high -> called Minus/high; margin 28.8 -> 28.8; best score 48.6 -> 48.6; models sexP:removed; genes algA,glrA,rnhA,sexM,tptA -> algA,glrA,sexM,tptA
- **GCA_977110945.1_gzUmbRama1** Mucoromycota:MAT CDSBDH010000020.1 [span_changed,gene_set_changed]: called Plus/high -> called Plus/high; margin 147.4 -> 147.4; best score 188.1 -> 188.1; genes algA,glrA,rnhA,sexP,tptA -> algA,glrA,sexP,tptA
- **GCA_977110975.1_gzUmbVina2** Mucoromycota:MAT CDSBDG010000016.1 [span_changed,gene_set_changed]: called Minus/high -> called Minus/high; margin 34.7 -> 34.7; best score 57.5 -> 57.5; genes algA,glrA,rnhA,sexM,sexP,tptA -> algA,glrA,sexM,tptA
- **GCF_025399195.1_Umbra1** Mucoromycota:MAT NW_026252087.1 [span_changed,gene_set_changed]: called Plus/high -> called Plus/high; margin 127.3 -> 127.3; best score 162.9 -> 162.9; genes algA,glrA,rnhA,sexP,tptA -> algA,glrA,sexP,tptA
