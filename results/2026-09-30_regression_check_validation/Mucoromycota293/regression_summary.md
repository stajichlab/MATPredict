# Regression check: curation-umbelopsis a2fe1b4 vs PR #9 8eff88e -- Mucoromycota293

- Genomes compared: 292; loci rows: 1820; loci with a change: 827 (134 touch a call, 693 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 656 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| call_gained | 1 |
| call_lost | 1 |
| classifier_input_changed | 9 |
| classifier_shift | 117 |
| confidence_changed | 8 |
| core_model_changed | 19 |
| gene_set_changed | 9 |
| idiomorph_changed | 5 |
| locus_class_changed | 8 |
| span_changed | 52 |

## Changes to calls

- **GCA_000534915.1_ASM53491v1** Mucoromycota:MAT BAVE01000001.1 [span_changed,core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 48.5 -> 78.7; best score 82.3 -> 119.8; models sexP:150bp/48.0%->783bp/31.13%
- **GCA_000587855.1_B50** Mucoromycota:MAT KK076501.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 30.5 -> 39.8; best score 81.8 -> 90.1
- **GCA_000696955.1_SynRacB6101-1.0** Mucoromycota:MAT JNDN01000979.1 [span_changed,core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 126.7 -> 189.7; best score 164.2 -> 227.7; models sexP:423bp/35.0%->930bp/100.0%
- **GCA_000696995.1_ApoEleB7760-1.0** Mucoromycota:MAT JNDQ01001420.1 [gene_set_changed]: called Minus/high -> called Minus/high; margin 84.4 -> 80.5; best score 167.8 -> 162.6; genes algA,sexM,tptA -> algA,sexM,sexP,tptA
- **GCA_000697035.1_RhiStoB9770-1.0** Mucoromycota:MAT JNDS01005206.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 274.7 -> 280.2; best score 329.6 -> 329.0
- **GCA_000697095.1_RhiDel18148-1.0** Mucoromycota:MAT JNDV01014111.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.1 -> 304.5; best score 359.0 -> 359.4
- **GCA_000697115.1_RhiOry21396-1.0** Mucoromycota:MAT JNDW01002981.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 298.3 -> 305.6; best score 358.5 -> 359.0
- **GCA_000697135.1_RhiOry99-133-1.0** Mucoromycota:MAT JNDX01004072.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_000697155.1_RhiDel21789-1.0** Mucoromycota:MAT JNDY01003786.1 [span_changed]: called Plus/high -> called Plus/high; margin 72.8 -> 74.2; best score 73.5 -> 74.2
- **GCA_000697195.1_MucRam97-1192-1.0** Mucoromycota:MAT JNEF01003481.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_000697235.1_CokRecB5483-1.0** Mucoromycota:MAT JNEH01001784.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 118.2 -> 112.2; best score 192.9 -> 187.2
- **GCA_000697255.1_MucRacB9645-1.0** Mucoromycota:MAT JNEI01002668.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 269.6 -> 262.4; best score 330.7 -> 327.9
- **GCA_000697275.1_RhiMicB9738-1.0** Mucoromycota:MAT JNEJ01004598.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 269.6 -> 262.4; best score 330.7 -> 327.9
- **GCA_000697355.1_SynMonB8922-1.0** Mucoromycota:MAT JNEN01001244.1 [span_changed,classifier_shift]: called Minus/medium -> called Minus/medium; margin 29.4 -> 35.1; best score 105.5 -> 108.9
- **GCA_000697415.1_UmbIsaB7317-1.0** Mucoromycota:MAT JNEQ01000056.1 [confidence_changed,locus_class_changed,span_changed,core_model_changed,classifier_input_changed]: called Plus/low -> called Plus/high; margin 92.3 -> 94.2; best score 121.9 -> 124.2; models sexP:0bp/42.857%->423bp/33.33%
- **GCA_000697435.1_RhiVarB7584-1.0** Mucoromycota:MAT JNES01000944.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 30.5 -> 39.8; best score 81.8 -> 90.1
- **GCA_000697605.1_RhiOryHUMC02-1.0** Mucoromycota:MAT KK960746.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 298.3 -> 305.6; best score 358.5 -> 359.0
- **GCA_000738595.1_RhiDel21447-1.0** Mucoromycota:MAT KL996174.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 291.2 -> 301.6; best score 354.8 -> 354.2
- **GCA_000738605.1_RhiDel21446-1.0** Mucoromycota:MAT KL995476.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 291.2 -> 301.6; best score 354.8 -> 354.2
- **GCA_001276145.1_ASM127614v1** Mucoromycota:MAT KQ435617.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 87.6 -> 76.7; best score 128.5 -> 117.0
- **GCA_001599575.1_JCM_22480_assembly_v001** Mucoromycota:MAT BCHG01000157.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 87.1 -> 74.4; best score 128.8 -> 116.9
- **GCA_001650995.1_ASM165099v1** Mucoromycota:MAT KV441911.1 [span_changed]: called Minus/high -> called Minus/high; margin 51.5 -> 49.4; best score 69.6 -> 67.1
- **GCA_002104935.1_Hesve2** Mucoromycota:MAT MCGT01000028.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 48.0 -> 57.2; best score 60.5 -> 65.8
- **GCA_002105135.1_Synrac1** Mucoromycota:MAT MCGN01000004.1 [span_changed,core_model_changed,classifier_shift]: called Plus/medium -> called Plus/medium; margin 126.7 -> 189.7; best score 164.2 -> 227.7; models sexP:423bp/35.0%->930bp/100.0%
- **GCA_002749535.1_ASM274953v1** Mucoromycota:MAT MZZL01000138.1 [gene_set_changed,classifier_shift]: called Minus/high -> called Minus/high; margin 85.7 -> 80.6; best score 169.4 -> 163.0; genes algA,sexM,tptA -> algA,sexM,sexP,tptA
- **GCA_003325435.1_Razy_CA** Mucoromycota:MAT PJQL01000316.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 147.0 -> 153.8; best score 216.6 -> 222.3
- **GCA_006680115.1_ASM668011v1** Mucoromycota:MAT SMRR01000092.1 [span_changed]: called Plus/high -> called Plus/high; margin 266.9 -> 263.2; best score 332.8 -> 329.7
- **GCA_010203745.1_Muccir1_3** Mucoromycota:MAT JAAECE010000001.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 125.4 -> 131.6; best score 171.5 -> 177.3
- **GCA_011763605.1_ASM1176360v1** Mucoromycota:MAT JAANIO010000501.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 274.7 -> 280.2; best score 329.6 -> 329.0
- **GCA_011763645.1_ASM1176364v1** Mucoromycota:MAT JAASGZ010000738.1 [span_changed]: called Minus/high -> called Minus/high; margin 174.2 -> 174.7; best score 245.4 -> 243.7
- **GCA_011763695.1_ASM1176369v1** Mucoromycota:MAT JAANIQ010000195.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 291.2 -> 301.6; best score 354.8 -> 354.2
- **GCA_011763715.1_ASM1176371v1** Mucoromycota:MAT JAANIR010000203.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 291.2 -> 301.6; best score 354.8 -> 354.2
- **GCA_011763755.1_ASM1176375v1** Mucoromycota:MAT JAANIW010000329.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 291.2 -> 301.6; best score 354.8 -> 354.2
- **GCA_011763835.1_ASM1176383v1** Mucoromycota:MAT JAANIX010000238.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 291.2 -> 301.6; best score 354.8 -> 354.2
- **GCA_011763935.1_ASM1176393v1** Mucoromycota:MAT JAANJC010000844.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 298.4 -> 306.2; best score 358.4 -> 359.0
- **GCA_011763955.1_ASM1176395v1** Mucoromycota:MAT JAANJD010000832.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 298.4 -> 306.2; best score 358.4 -> 359.0
- **GCA_011763975.1_ASM1176397v1** Mucoromycota:MAT JAANJF010000979.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 298.7 -> 305.5; best score 359.8 -> 360.3
- **GCA_011764055.1_ASM1176405v1** Mucoromycota:MAT JAANQS010002783.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011764145.1_ASM1176414v1** Mucoromycota:MAT JAANQZ010003255.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.3 -> 304.9; best score 359.5 -> 359.9
- **GCA_011764185.1_ASM1176418v1** Mucoromycota:MAT JAANRA010001051.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011764225.1_ASM1176422v1** Mucoromycota:MAT JAANRQ010003168.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011800955.2_ASM1180095v2** Mucoromycota:MAT JAANRM010001105.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011800985.2_ASM1180098v2** Mucoromycota:MAT JAANRL010002477.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011801035.2_ASM1180103v2** Mucoromycota:MAT JAANRN010000391.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011801055.2_ASM1180105v2** Mucoromycota:MAT JAANRO010000810.1 [span_changed,classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011801305.2_ASM1180130v2** Mucoromycota:MAT JAANRS010001000.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011801385.2_ASM1180138v2** Mucoromycota:MAT JAANRR010002648.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011801495.1_ASM1180149v1** Mucoromycota:MAT JAANRV010002415.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011801505.2_ASM1180150v2** Mucoromycota:MAT JAANRU010002094.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011801555.2_ASM1180155v2** Mucoromycota:MAT JAANRX010000827.1 [span_changed,classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011801575.2_ASM1180157v2** Mucoromycota:MAT JAANRW010000843.1 [span_changed,classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011801605.2_ASM1180160v2** Mucoromycota:MAT JAANRT010001072.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011801645.2_ASM1180164v2** Mucoromycota:MAT JAANRB010001607.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011801705.2_ASM1180170v2** Mucoromycota:MAT JAANQY010002666.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011801725.2_ASM1180172v2** Mucoromycota:MAT JAANQX010000106.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011801815.2_ASM1180181v2** Mucoromycota:MAT JAANQU010000727.1 [span_changed,classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011950945.2_ASM1195094v2** Mucoromycota:MAT JAANRD010001073.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011950985.2_ASM1195098v2** Mucoromycota:MAT JAANRF010000404.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 297.6 -> 304.4; best score 358.8 -> 359.2
- **GCA_011951295.2_ASM1195129v2** Mucoromycota:MAT JAANRC010000468.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011952115.2_ASM1195211v2** Mucoromycota:MAT JAANRI010000651.1 [span_changed,classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011952145.2_ASM1195214v2** Mucoromycota:MAT JAANRJ010002406.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_011952175.2_ASM1195217v2** Mucoromycota:MAT JAANRK010002299.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 298.7 -> 305.3; best score 359.4 -> 359.9
- **GCA_013461545.1_ASM1346154v1** Mucoromycota:MAT VAFG01000366.1 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 29.7 -> 35.1; best score 102.0 -> 105.7
- **GCA_016758895.1_ASM1675889v1** Mucoromycota:MAT JAEPRA010000016.1 [idiomorph_changed,confidence_changed,locus_class_changed,span_changed,gene_set_changed,core_model_changed,classifier_input_changed,classifier_shift]: called undetermined/low -> called Minus/high; margin 13.0 -> 100.6; best score 41.0 -> 127.3; models sexM:0bp/36.232%->450bp/100.0%; genes algA,glrA,sexM,tptA -> algA,glrA,rnhA,sexM,tptA
- **GCA_016758905.1_ASM1675890v1** Mucoromycota:MAT JAEPQZ010000017.1 [confidence_changed,locus_class_changed,span_changed,gene_set_changed,core_model_changed,classifier_input_changed,classifier_shift]: called Plus/low -> called Plus/medium; margin 68.1 -> 76.5; best score 100.4 -> 111.2; models sexP:0bp/41.667%->462bp/32.45%; genes algA,sexP,tptA -> algA,sexM,sexP,tptA
- **GCA_016758965.1_ASM1675896v1** Mucoromycota:MAT JAEPRB010000020.1 [call_lost,classifier_shift]: called Plus/medium -> withheld:mat_gene_gate Plus/; margin 81.7 -> 54.2; best score 111.8 -> 86.4; models no_gene_evidence_in_candidate
- **GCA_019022795.1_ASM1902279v1** Mucoromycota:MAT JAFCNC010000323.1 [span_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 291.2 -> 301.6; best score 354.8 -> 354.2
- **GCA_019677225.1_ASM1967722v1** Mucoromycota:MAT JAFHAR010000822.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 58.8 -> 51.1; best score 112.3 -> 102.6
- **GCA_021235465.1_Umbelo1** Mucoromycota:MAT JAJSME010000005.1 [core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 61.6 -> 135.3; best score 86.8 -> 165.6; models sexP:237bp/40.79%->915bp/56.77%
- **GCA_022702545.1_SJTU-YXQ** Mucoromycota:MAT JAHRBC010000002.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 131.5 -> 125.8; best score 212.6 -> 206.1
- **GCA_023629755.1_ASM2362975v1** Mucoromycota:MAT JAMAMJ010000480.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 87.6 -> 76.7; best score 128.5 -> 117.0
- **GCA_023629815.1_ASM2362981v1** Mucoromycota:MAT JAMAMK010000501.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 87.6 -> 76.7; best score 128.5 -> 117.0
- **GCA_023630295.1_ASM2363029v1** Mucoromycota:MAT JAMAMR010000041.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 298.1 -> 304.5; best score 359.0 -> 359.4
- **GCA_023630305.1_ASM2363030v1** Mucoromycota:MAT JAMAMS010000120.1 [span_changed,core_model_changed,classifier_shift]: called Plus/medium -> called Plus/medium; margin 131.1 -> 164.2; best score 165.3 -> 195.8; models sexP:501bp/31.93%->918bp/75.16%
- **GCA_023630315.1_ASM2363031v1** Mucoromycota:MAT JAMAMU010000151.1 [span_changed,core_model_changed,classifier_shift]: called Plus/medium -> called Plus/medium; margin 122.0 -> 158.4; best score 158.0 -> 190.6; models sexP:420bp/35.25%->933bp/72.58%
- **GCA_023630325.1_ASM2363032v1** Mucoromycota:MAT JAMAMT010000060.1 [span_changed,core_model_changed,classifier_shift]: called Plus/medium -> called Plus/medium; margin 122.0 -> 158.4; best score 158.0 -> 190.6; models sexP:420bp/35.25%->933bp/72.58%
- **GCA_023630435.1_ASM2363043v1** Mucoromycota:MAT JAMAMP010001032.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 84.1 -> 71.9; best score 123.8 -> 110.4
- **GCA_023630465.1_ASM2363046v1** Mucoromycota:MAT JAMAMO010000771.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 84.1 -> 71.9; best score 123.8 -> 110.4
- **GCA_023630475.1_ASM2363047v1** Mucoromycota:MAT JAMAMN010001208.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 84.1 -> 71.9; best score 123.8 -> 110.4
- **GCA_023630555.1_ASM2363055v1** Mucoromycota:MAT JAMAMM010000510.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 87.6 -> 76.7; best score 128.5 -> 117.0
- **GCA_024139275.1_ASM2413927v1** Mucoromycota:MAT VCJA01002019.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 243.0 -> 252.7; best score 310.7 -> 318.4
- **GCA_024139395.1_ASM2413939v1** Mucoromycota:MAT VCJI01001366.1 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 75.8 -> 94.0; best score 150.6 -> 166.7
- **GCA_024139405.1_ASM2413940v1** Mucoromycota:MAT VCJL01001223.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 77.0 -> 89.3; best score 153.5 -> 163.8
- **GCA_024140155.1_ASM2414015v1** Mucoromycota:MAT VCJB01002633.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 256.5 -> 263.0; best score 321.9 -> 326.5
- **GCA_024140195.1_ASM2414019v1** Mucoromycota:MAT VCJC01002502.1 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 74.4 -> 92.0; best score 150.1 -> 166.3
- **GCA_024140245.1_ASM2414024v1** Mucoromycota:MAT VCJD01005173.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 255.9 -> 265.4; best score 324.0 -> 331.5
- **GCA_024140275.1_ASM2414027v1** Mucoromycota:MAT VCJF01002266.1 [span_changed,classifier_shift]: called Minus/high -> called Minus/high; margin 76.8 -> 92.2; best score 152.1 -> 166.8
- **GCA_024140565.1_ASM2414056v1** Mucoromycota:MAT VCJH01002494.1 [span_changed]: called Plus/high -> called Plus/high; margin 221.4 -> 220.7; best score 292.9 -> 294.7
- **GCA_025201375.1_Gonbut1** Mucoromycota:MAT JAIWZP010000006.1 [span_changed]: called Minus/high -> called Minus/high; margin 50.7 -> 48.0; best score 72.1 -> 69.1
- **GCA_025201805.1_Parpar1** Mucoromycota:MAT MU619237.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 130.1 -> 123.9; best score 174.5 -> 166.9
- **GCA_025266875.1_Fenlin1** Mucoromycota:MAT JALLLT010000045.1 [span_changed]: called Minus/medium -> called Minus/medium; margin 56.7 -> 56.4; best score 122.9 -> 121.8
- **GCA_025331445.1_Blatri1** Mucoromycota:MAT JAIWNE010000009.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 239.8 -> 246.2; best score 301.7 -> 299.8
- **GCA_025399235.1_Bacci1** Mucoromycota:MAT MU622345.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 239.3 -> 248.3; best score 304.6 -> 312.8
- **GCA_025504515.1_ASM2550451v1** Mucoromycota:MAT JAMFRH010000492.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 87.6 -> 76.7; best score 128.5 -> 117.0
- **GCA_025529335.1_ASM2552933v1** Mucoromycota:MAT JAMYHX010000536.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 120.5 -> 125.5; best score 193.6 -> 197.3
- **GCA_025529335.1_ASM2552933v1** Mucoromycota:MAT JAMYHX010000515.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 272.2 -> 265.9; best score 327.5 -> 324.6
- **GCA_027595865.1_Umbelopsis_isabellina_MPG-14A_reass_enz** Mucoromycota:MAT JANJFK010000028.1 [confidence_changed,locus_class_changed,span_changed,core_model_changed,classifier_input_changed]: called Plus/low -> called Plus/high; margin 92.3 -> 94.2; best score 121.9 -> 124.2; models sexP:0bp/42.857%->423bp/33.33%
- **GCA_037039815.1_ASM3703981v1** Mucoromycota:MAT JAXQGL010000038.1 [confidence_changed,locus_class_changed,span_changed,gene_set_changed,core_model_changed,classifier_input_changed]: called Plus/low -> called Plus/high; margin 75.8 -> 77.8; best score 103.7 -> 106.1; models sexP:0bp/33.333%->777bp/28.91%; genes algA,sexP,tptA -> algA,glrA,sexM,sexP,tptA
- **GCA_037042155.1_ASM3704215v1** Mucoromycota:MAT JAXQJX010000258.1 [span_changed,core_model_changed,classifier_shift]: called Plus/medium -> called Plus/medium; margin 131.1 -> 164.3; best score 165.3 -> 195.6; models sexP:501bp/31.93%->918bp/75.16%
- **GCA_039881115.1_ASM3988111v1** Mucoromycota:MAT JARGYK010002198.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 242.0 -> 249.3; best score 300.8 -> 303.7
- **GCA_040206765.1_ASM4020676v1** Mucoromycota:MAT JBEEFM010000005.1 [span_changed,classifier_shift]: called Minus/high -> called Minus/high; margin 120.5 -> 125.5; best score 193.6 -> 197.3
- **GCA_040256725.1_ASM4025672v1** Mucoromycota:MAT MU975357.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 124.3 -> 132.7; best score 171.0 -> 178.9
- **GCA_049863895.1_ASM4986389v1** Mucoromycota:MAT SMHB01000855.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 299.1 -> 305.9; best score 359.1 -> 359.5
- **GCA_050656845.1_ASM5065684v1** Mucoromycota:MAT JBOEPV010000879.1 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 62.3 -> 57.3; best score 106.9 -> 100.0
- **GCA_051997045.1_ASM5199704v1** Mucoromycota:MAT JARGYJ010000412.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 145.8 -> 151.7; best score 215.7 -> 220.3
- **GCA_051997065.1_ASM5199706v1** Mucoromycota:MAT JARGYI010000307.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 145.8 -> 151.7; best score 215.7 -> 220.3
- **GCA_051997145.1_ASM5199714v1** Mucoromycota:MAT JARGYL010000220.1 [span_changed]: called Minus/high -> called Minus/high; margin 138.5 -> 140.9; best score 210.9 -> 211.7
- **GCA_051997165.1_ASM5199716v1** Mucoromycota:MAT JARGYM010000377.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 145.8 -> 151.7; best score 215.7 -> 220.3
- **GCA_053572175.1_ASM5357217v1** Mucoromycota:MAT JAYKKX010000013.1 [idiomorph_changed,confidence_changed,locus_class_changed,span_changed,core_model_changed,classifier_input_changed,classifier_shift]: called undetermined/low -> called Minus/high; margin 16.3 -> 37.1; best score 33.8 -> 48.3; models sexM:0bp/36.232%->189bp/47.62%
- **GCA_054601915.1_FAM010_1.0** Mucoromycota:MAT BAAHPX010000008.1 [span_changed]: called Minus/high -> called Minus/high; margin 77.3 -> 76.2; best score 133.0 -> 130.4
- **GCA_054906655.1_Umbelopsis_nana_v._1.0** Mucoromycota:MAT BAAHSD010000001.1 [call_gained,gene_set_changed,classifier_input_changed,classifier_shift]: withheld:flank_carried_core_outside_flank_span undetermined/ -> called Minus/high; margin 6.4 -> 31.2; best score 25.2 -> 50.8; models no_gene_evidence_in_baseline; genes algA,glrA,sexP,tptA -> algA,glrA,sexM,tptA
- **GCA_054954275.1_ASM5495427v1** Mucoromycota:MAT CM143249.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 109.4 -> 117.9; best score 176.7 -> 183.3
- **GCA_057662165.1_ASM5766216v1** Mucoromycota:MAT JBYEQD010000037.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 289.5 -> 301.1; best score 354.1 -> 353.3
- **GCA_058380205.1_DTU-UU_Roli_1.0** Mucoromycota:MAT JBYHKW010000001.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 269.6 -> 262.4; best score 330.7 -> 327.9
- **GCA_059058645.1_ASM5905864v1** Mucoromycota:MAT CM180848.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 109.4 -> 117.9; best score 176.7 -> 183.3
- **GCA_059714295.1_ASM5971429v1** Mucoromycota:MAT JCAOFY010000002.1 [span_changed,core_model_changed,classifier_shift]: called Plus/medium -> called Plus/medium; margin 119.0 -> 154.3; best score 155.1 -> 190.7; models sexP:417bp/33.09%->927bp/78.1%
- **GCA_060230415.1_ASM6023041v1** Mucoromycota:MAT CM190389.1 [span_changed,classifier_shift]: called Minus/high -> called Minus/high; margin 120.5 -> 125.5; best score 193.6 -> 197.3
- **GCA_060230415.1_ASM6023041v1** Mucoromycota:MAT CM190391.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 269.6 -> 262.4; best score 330.7 -> 327.9
- **GCA_060309175.1_CBS116.08v1** Mucoromycota:MAT JBLNHG010000007.1 [idiomorph_changed,classifier_shift]: called Minus/high -> called undetermined/high; margin 26.3 -> 17.8; best score 31.7 -> 25.0
- **GCA_060309305.1_CBS205.68v1** Mucoromycota:MAT JBLNHI010000055.1 [span_changed]: called Plus/high -> called Plus/high; margin 319.4 -> 320.1; best score 372.8 -> 369.2
- **GCA_060309335.1_CBS210.80v1** Mucoromycota:MAT JBLNHJ010000021.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 86.6 -> 93.3; best score 158.4 -> 162.1
- **GCA_900079185.1_AG_v1** Mucoromycota:MAT LT553527.1 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 107.5 -> 91.1; best score 158.2 -> 141.7
- **GCA_964291815.1_UMBE_WA70503** Mucoromycota:MAT CAXYTU010000169.1 [idiomorph_changed,confidence_changed,locus_class_changed,span_changed,core_model_changed,classifier_input_changed,classifier_shift]: called undetermined/low -> called Minus/high; margin 16.3 -> 37.1; best score 33.8 -> 48.3; models sexM:0bp/36.232%->189bp/47.62%
- **GCA_977092135.1_gzAbsCyli1** Mucoromycota:MAT CDRYNH010000023.1 [classifier_shift]: called Minus/medium -> called Minus/medium; margin 51.6 -> 46.1; best score 111.5 -> 106.0
- **GCA_977110905.1_gzMucHiem1** Mucoromycota:MAT CDSBDF010000023.1 [classifier_shift]: called Minus/low -> called Minus/low; margin 45.5 -> 52.1; best score 97.7 -> 103.7
- **GCA_977110945.1_gzUmbRama1** Mucoromycota:MAT CDSBDH010000020.1 [span_changed,gene_set_changed,core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 56.0 -> 137.9; best score 87.8 -> 169.7; models sexP:240bp/30.0%->960bp/100.0%; genes algA,glrA,sexP,tptA -> algA,glrA,rnhA,sexP,tptA
- **GCA_977110975.1_gzUmbVina2** Mucoromycota:MAT CDSBDG010000016.1 [idiomorph_changed,confidence_changed,locus_class_changed,span_changed,gene_set_changed,core_model_changed,classifier_input_changed,classifier_shift]: called undetermined/low -> called Minus/high; margin 8.9 -> 36.8; best score 31.8 -> 58.3; models sexM:0bp/34.783%->378bp/53.97%; genes algA,sexM,tptA -> algA,glrA,sexM,sexP,tptA
- **GCF_000149305.1_RO3** Mucoromycota:MAT NW_027137827.1 [classifier_shift]: called Plus/high -> called Plus/high; margin 291.2 -> 301.6; best score 354.8 -> 354.2
- **GCF_025118155.1_Cokrec1** Mucoromycota:MAT NW_026251389.1 [classifier_shift]: called Minus/high -> called Minus/high; margin 118.2 -> 112.2; best score 192.9 -> 187.2
- **GCF_025118155.1_Cokrec1** Mucoromycota:MAT NW_026251417.1 [span_changed]: called Plus/medium -> called Plus/medium; margin 294.4 -> 293.9; best score 353.6 -> 352.5
- **GCF_025201335.1_Gilper1** Mucoromycota:MAT NW_026251964.1 [span_changed]: called Minus/high -> called Minus/high; margin 77.3 -> 76.2; best score 133.0 -> 130.4
- **GCF_025201355.1_Halrad1** Mucoromycota:MAT NW_026251697.1 [span_changed]: called Minus/high -> called Minus/high; margin 49.2 -> 54.0; best score 110.2 -> 111.0
- **GCF_025331425.1_Radspe1** Mucoromycota:MAT NW_026251940.1 [classifier_shift]: called Plus/medium -> called Plus/medium; margin 179.9 -> 172.9; best score 314.8 -> 307.7
- **GCF_025399195.1_Umbra1** Mucoromycota:MAT NW_026252087.1 [span_changed,gene_set_changed,core_model_changed,classifier_shift]: called Plus/high -> called Plus/high; margin 51.1 -> 126.7; best score 81.6 -> 158.3; models sexP:264bp/29.55%->951bp/84.71%; genes algA,glrA,sexP,tptA -> algA,glrA,rnhA,sexP,tptA
