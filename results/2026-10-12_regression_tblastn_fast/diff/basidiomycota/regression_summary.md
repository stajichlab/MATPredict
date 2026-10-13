# Regression check: regression run-9b15d43 vs run-73f56b8 -- basidiomycota

- Genomes compared: 33; loci rows: 4055; loci with a change: 3533 (16 touch a call, 3517 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 0 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| call_lost | 2 |
| core_model_changed | 4 |
| gene_set_changed | 4 |
| span_changed | 14 |

## Changes to calls

- **GCA_000715385.1_Rhizoctonia_solani_123E_v1** Basidiomycota:PR AZST01000319.1 [span_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> 
- **GCA_001542265.1_ASM154226v1** Basidiomycota:redHD LNKU01000002.1 [span_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> 
- **GCA_002092955.1_ASM209295v1** Basidiomycota:PR NCVV01000007.1 [span_changed,gene_set_changed,core_model_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> ; models fungal_mating_type_pheromone:removed; genes caax_precursor,fungal_mating_type_pheromone,pheromone_receptor -> caax_precursor,pheromone_receptor
- **GCA_019143615.1_ASM1914361v1** Basidiomycota:Aalpha JAGVSI010001513.1 [span_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> 
- **GCA_019143615.1_ASM1914361v1** Basidiomycota:HD JAGVSI010001121.1 [span_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> 
- **GCA_019143615.1_ASM1914361v1** Basidiomycota:Balpha JAGVSI010000976.1 [span_changed,gene_set_changed,core_model_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; models bbp2-1:removed;bbp2-8:removed;fungal_mating_type_pheromone:removed; genes bap3-1,bap3-3,bar3,bbp2-1,bbp2-2,bbp2-5,bbp2-6,bbp2-8,bbp2_b,bbr2,caax_precursor,fungal_mating_type_pheromone,pheromone_receptor -> bap3-1,bap3-3,bar3,bbp2-2,bbp2-5,bbp2-6,bbp2_b,bbr2,caax_precursor,pheromone_receptor
- **GCA_023273805.1_ASM2327380v1** Basidiomycota:MAT CP096880.1 [call_lost,span_changed]: called undetermined/medium -> withheld:below_fraction_floor undetermined/; margin  -> ; best score  -> ; models no_gene_evidence_in_candidate
- **GCA_023273805.1_ASM2327380v1** Basidiomycota:bLocus CP096878.1 [span_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> 
- **GCA_050403325.1_S_johnsonii_JCM1840v1** Basidiomycota:redHD BAABVI010000025.1 [span_changed,core_model_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; models HD1:876bp/42.91%->867bp/43.71%
- **GCA_056247695.1_ASM5624769v1** Basidiomycota:PR JBVYWQ010001089.1 [call_lost,gene_set_changed]: called undetermined/medium -> withheld:below_fraction_floor undetermined/; margin  -> ; best score  -> ; models no_gene_evidence_in_candidate; genes caax_precursor,fungal_mating_type_pheromone,pheromone_receptor -> caax_precursor,pheromone_receptor
- **GCA_900162835.1_version_2** Basidiomycota:aLocus LT795073.1 [span_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> 
- **GCA_943193715.1_gfAgaBisp1.1** Basidiomycota:HD OW971906.1 [span_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> 
- **GCF_000264905.1_Stehi1** Basidiomycota:HD NW_006763131.1 [span_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> 
- **GCF_000271585.1_Trametes_versicolor_v1.0** Basidiomycota:HD NW_007360321.1 [span_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> 
- **GCF_000271585.1_Trametes_versicolor_v1.0** Basidiomycota:PR NW_007360328.1 [gene_set_changed,core_model_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> ; models pheromone_B43:removed; genes caax_precursor,fungal_mating_type_pheromone,pheromone_B43,pheromone_receptor -> caax_precursor,fungal_mating_type_pheromone,pheromone_receptor
- **GCF_000271605.1_Fomme1** Basidiomycota:HD NW_006760390.1 [span_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> 
