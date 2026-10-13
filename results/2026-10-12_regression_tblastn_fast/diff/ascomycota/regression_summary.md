# Regression check: regression run-9b15d43 vs run-73f56b8 -- ascomycota

- Genomes compared: 30; loci rows: 1181; loci with a change: 1133 (16 touch a call, 1117 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 0 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| confidence_changed | 2 |
| core_model_changed | 9 |
| gene_set_changed | 7 |
| span_changed | 11 |

## Changes to calls

- **GCA_000147375.2_WaNRRLY366_v02** Ascomycota:MATsc JH119216.1 [span_changed]: called MATa/medium -> called MATa/medium; margin  -> ; best score  -> 
- **GCA_000147375.2_WaNRRLY366_v02** Ascomycota:MATsc JH119233.1 [gene_set_changed,core_model_changed]: called MATalpha/medium -> called MATalpha/medium; margin  -> ; best score  -> ; models MATA2:207bp/41.54%->0bp/40.625%; genes MATA2,MATALPHA1,MATALPHA2,sla2 -> MATALPHA1,MATALPHA2,sla2
- **GCA_001510565.1_ASM151056v1** Ascomycota:MAT LNGL01000024.1 [span_changed]: called MAT1-2/medium -> called MAT1-2/medium; margin  -> ; best score  -> 
- **GCA_002924885.1_ASM292488v1** Ascomycota:MAT PJYB01000136.1 [span_changed]: called MAT1-1/high -> called MAT1-1/high; margin  -> ; best score  -> 
- **GCA_003705895.2_ASM370589v2** Ascomycota:MAT PPJV02000047.1 [confidence_changed,span_changed,gene_set_changed,core_model_changed]: called MAT1-2/high -> called MAT1-2/medium; margin  -> ; best score  -> ; models MAT1-1-1:removed;MAT1-1-2:removed;MAT1-1-4:removed; genes APN2,MAT1-1-1,MAT1-1-2,MAT1-1-4,MAT1-2-1,SLA2 -> APN2,MAT1-2-1,SLA2
- **GCA_003705895.2_ASM370589v2** Ascomycota:MATtub PPJV02000047.1 [span_changed]: called MAT1-2/medium -> called MAT1-2/medium; margin  -> ; best score  -> 
- **GCA_004918105.1_ASM491810v1** Ascomycota:MAT QZBS01000102.1 [span_changed,gene_set_changed]: called MAT1-2/medium -> called MAT1-2/medium; margin  -> ; best score  -> ; genes APN2,COX13,MAT1-1-1,MAT1-2-1,SLA2 -> APN2,COX13,MAT1-1-1,MAT1-2-1
- **GCA_010015785.1_Sacpr1** Ascomycota:MAT ML978713.1 [confidence_changed,gene_set_changed,core_model_changed]: called MAT1-2/medium -> called MAT1-2/high; margin  -> ; best score  -> ; models MAT1-1-1:removed;MAT1-1-4:removed; genes APN2,COX13,MAT1-1-4,MAT1-2-1 -> APN2,COX13,MAT1-2-1
- **GCA_022813925.1_ASM2281392v1** Ascomycota:MAT JALDZU010000168.1 [span_changed,gene_set_changed,core_model_changed]: called MAT1-2/high -> called MAT1-2/high; margin  -> ; best score  -> ; models MAT1-1-1:removed;MAT1-1-3:removed; genes APN2,COX13,MAT1-1-1,MAT1-1-3,MAT1-2-1 -> APN2,COX13,MAT1-2-1
- **GCA_022814105.1_ASM2281410v1** Ascomycota:MAT JALDZM010000006.1 [span_changed,gene_set_changed,core_model_changed]: called MAT1-2/high -> called MAT1-2/high; margin  -> ; best score  -> ; models MAT1-1-1:removed; genes APN2,COX13,MAT1-1-1,MAT1-2-1,SLA2 -> APN2,COX13,MAT1-2-1,SLA2
- **GCA_023645945.1_KEW_Ddia** Ascomycota:MAT JALAIB010000011.1 [span_changed]: called MAT1-1/medium -> called MAT1-1/medium; margin  -> ; best score  -> 
- **GCA_025802115.2_ASM2580211v2** Ascomycota:MAT JAOYNC020110414.1 [span_changed,gene_set_changed,core_model_changed]: called MAT1-1/medium -> called MAT1-1/medium; margin  -> ; best score  -> ; models MAT1-1-2:removed; genes MAT1-1-1,MAT1-1-2,SLA2 -> MAT1-1-1,SLA2
- **GCA_026023515.1_ASM2602351v1** Ascomycota:MAT JAPETT010000019.1 [span_changed]: called MAT1-1/medium -> called MAT1-1/medium; margin  -> ; best score  -> 
- **GCA_026023555.1_ASM2602355v1** Ascomycota:MAT JAPETR010000115.1 [core_model_changed]: called MAT1-2/high -> called MAT1-2/high; margin  -> ; best score  -> ; models MAT1-1-1:removed
- **GCA_964254465.1_private_X5_concoct_bin.21** Ascomycota:MAT CAXUUA010000001.1 [core_model_changed]: called MAT1-2/high -> called MAT1-2/high; margin  -> ; best score  -> ; models MAT1-1-1:removed
- **GCA_964255695.1_public_SRR14722165_metabat2_bin.16** Ascomycota:MAT CAXUYV010000011.1 [core_model_changed]: called MAT1-2/high -> called MAT1-2/high; margin  -> ; best score  -> ; models MAT1-1-1:removed
