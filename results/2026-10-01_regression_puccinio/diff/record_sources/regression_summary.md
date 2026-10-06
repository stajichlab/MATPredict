# Regression check: regression run-cand-puccinio vs run-c26669c -- record_sources

- Genomes compared: 23; loci rows: 911; loci with a change: 150 (4 touch a call, 146 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 0 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| core_model_changed | 1 |
| gene_set_changed | 3 |
| locus_class_changed | 1 |
| span_changed | 4 |

## Changes to calls

- **GCA_016772295.1_ASM1677229v1** Basidiomycota:Aalpha JAAGWA010000001.1 [span_changed,gene_set_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; genes HD1,HD2,Y,Z -> HD1,HD2,MIP1,Y,Z,beta_fg
- **GCF_000143185.2_Schco3** Basidiomycota:HD NW_026089539.1 [span_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> 
- **GCF_000143185.2_Schco3** Basidiomycota:Aalpha NW_026089539.1 [span_changed,gene_set_changed,core_model_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; models HD1:105bp/42.86%->255bp/34.52%; genes HD1,HD2,Y,Z -> HD1,HD2,MIP1,Y,Z
- **GCF_000300575.1_Agabi_varbisH97_2** Basidiomycota:HD NW_006267344.1 [locus_class_changed,span_changed,gene_set_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; genes HD1,HD2 -> HD1,HD2,MIP1,beta_fg
