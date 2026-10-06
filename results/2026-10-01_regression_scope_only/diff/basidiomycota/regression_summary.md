# Regression check: regression cand-fix vs run-cand-basidio -- basidiomycota

- Genomes compared: 33; loci rows: 5117; loci with a change: 1677 (19 touch a call, 1658 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 0 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| call_gained | 2 |
| call_lost | 15 |
| span_changed | 2 |

## Changes to calls

- **GCA_000715385.1_Rhizoctonia_solani_123E_v1** Basidiomycota:redPR AZST01000156.1 [call_lost]: called A2/high -> absent /; margin  -> ; best score  -> 
- **GCA_000715385.1_Rhizoctonia_solani_123E_v1** Basidiomycota:wallMAT AZST01000156.1 [call_lost]: called v1/medium -> absent /; margin  -> ; best score  -> 
- **GCA_000715385.1_Rhizoctonia_solani_123E_v1** Basidiomycota:wallMAT AZST01000319.1 [call_lost]: called v2/medium -> absent /; margin  -> ; best score  -> 
- **GCA_000715385.1_Rhizoctonia_solani_123E_v1** Basidiomycota:PR AZST01000319.1 [span_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> 
- **GCA_000715385.1_Rhizoctonia_solani_123E_v1** Basidiomycota:wallMAT AZST01000540.1 [call_lost]: called v1/medium -> absent /; margin  -> ; best score  -> 
- **GCA_000715385.1_Rhizoctonia_solani_123E_v1** Basidiomycota:PR AZST01000540.1 [span_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> 
- **GCA_001244265.2_MvSlA1A2r3c** Basidiomycota:redPR LN717261.1 [call_lost]: called A1/medium -> absent /; margin  -> ; best score  -> 
- **GCA_001244265.2_MvSlA1A2r3c** Basidiomycota:wallMAT LN717261.1 [call_lost]: called v2/medium -> absent /; margin  -> ; best score  -> 
- **GCA_002092955.1_ASM209295v1** Basidiomycota:redPR NCVV01000007.1 [call_lost]: called A1/high -> absent /; margin  -> ; best score  -> 
- **GCA_002092955.1_ASM209295v1** Basidiomycota:wallMAT NCVV01000007.1 [call_lost]: called v1/medium -> absent /; margin  -> ; best score  -> 
- **GCA_002105055.1_Leucr1** Basidiomycota:redHD MCGR01000028.1 [call_lost]: called undetermined/high -> absent /; margin  -> ; best score  -> 
- **GCA_002105055.1_Leucr1** Basidiomycota:redPR MCGR01000123.1 [call_lost]: called A2/medium -> absent /; margin  -> ; best score  -> 
- **GCA_023273805.1_ASM2327380v1** Basidiomycota:wallMAT CP096879.1 [call_lost]: called v1/medium -> absent /; margin  -> ; best score  -> 
- **GCF_000218685.1_v1.0** Basidiomycota:PR NW_006763292.1 [call_gained]: absent / -> called undetermined/medium; margin  -> ; best score  -> 
- **GCF_000218685.1_v1.0** Basidiomycota:PR NW_006763293.1 [call_gained]: absent / -> called undetermined/medium; margin  -> ; best score  -> 
- **GCF_000271605.1_Fomme1** Basidiomycota:redPR NW_006760389.1 [call_lost]: called A1/high -> absent /; margin  -> ; best score  -> 
- **GCF_000271605.1_Fomme1** Basidiomycota:wallMAT NW_006760389.1 [call_lost]: called v2/medium -> absent /; margin  -> ; best score  -> 
- **GCF_000292625.1_Dacryopinax_sp._DJM_731_SSP1_v1.0** Basidiomycota:redPR NW_024467210.1 [call_lost]: called A1/high -> absent /; margin  -> ; best score  -> 
- **GCF_000292625.1_Dacryopinax_sp._DJM_731_SSP1_v1.0** Basidiomycota:wallMAT NW_024467210.1 [call_lost]: called v2/medium -> absent /; margin  -> ; best score  -> 
