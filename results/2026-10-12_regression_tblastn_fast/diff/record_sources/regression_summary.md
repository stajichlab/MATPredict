# Regression check: regression run-9b15d43 vs run-73f56b8 -- record_sources

- Genomes compared: 23; loci rows: 966; loci with a change: 829 (10 touch a call, 819 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 0 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| confidence_changed | 1 |
| core_model_changed | 6 |
| gene_set_changed | 3 |
| locus_class_changed | 3 |
| span_changed | 6 |

## Changes to calls

- **GCA_023646395.1_KEW_Xsp1** Ascomycota:MAT JALAIM010000088.1 [locus_class_changed,core_model_changed]: called MAT1-2/high -> called MAT1-2/high; margin  -> ; best score  -> ; models MAT1-1-3:removed
- **GCA_023647215.1_KEW_Xsp2** Ascomycota:MAT JALAIR010000062.1 [locus_class_changed,core_model_changed]: called MAT1-2/high -> called MAT1-2/high; margin  -> ; best score  -> ; models MAT1-1-3:removed
- **GCA_964255735.1_public_SRR14722176_metabat2_bin.14** Ascomycota:MAT CAXUYH010000110.1 [core_model_changed]: called MAT1-2/medium -> called MAT1-2/medium; margin  -> ; best score  -> ; models MAT1-1-3:removed
- **GCA_964255735.1_public_SRR14722176_metabat2_bin.14** Ascomycota:MAT CAXUYH010000065.1 [span_changed,gene_set_changed,core_model_changed]: called MAT1-1/high -> called MAT1-1/high; margin  -> ; best score  -> ; models MAT1-1-2:removed;MAT1-1-3:removed; genes APN2,COX13,MAT1-1-1,MAT1-1-2,MAT1-1-3 -> APN2,COX13,MAT1-1-1
- **GCF_000002515.2_ASM251v1** Ascomycota:MATsc NC_006039.1 [span_changed]: called MATalpha/high -> called MATalpha/high; margin  -> ; best score  -> 
- **GCF_000002945.2_ASM294v3** Ascomycota:mat1 NC_003423.3 [span_changed]: called M/high -> called M/high; margin  -> ; best score  -> 
- **GCF_000142805.1_ASM14280v1** Ascomycota:MATsc NC_013082.1 [gene_set_changed,core_model_changed]: called MATalpha/high -> called MATalpha/high; margin  -> ; best score  -> ; models MATA2:210bp/40.0%->0bp/45.763%; genes MATA1,MATA2,MATALPHA1,MATALPHA2 -> MATA1,MATALPHA1,MATALPHA2
- **GCF_000142805.1_ASM14280v1** Ascomycota:MATsc NC_013082.1 [confidence_changed,locus_class_changed,span_changed,gene_set_changed,core_model_changed]: called MATa/high -> called MATa/medium; margin  -> ; best score  -> ; models MATA2:removed; genes MATA1,MATA2,sla2 -> MATA1,sla2
- **GCF_000143185.2_Schco3** Basidiomycota:HD NW_026089539.1 [span_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> 
- **GCF_000300575.1_Agabi_varbisH97_2** Basidiomycota:HD NW_006267344.1 [span_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> 
