# Regression check: net: PR#9 + puccinio + basidio-anchors + scope-only fix vs PR#9 -- basidiomycota

- Genomes compared: 33; loci rows: 6842; loci with a change: 4684 (35 touch a call, 4649 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 0 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| call_gained | 13 |
| call_lost | 2 |
| core_model_changed | 16 |
| gene_set_changed | 22 |
| locus_class_changed | 6 |
| span_changed | 24 |
| verification_changed | 3 |

## Changes to calls

- **GCA_000715385.1_Rhizoctonia_solani_123E_v1** Basidiomycota:PR AZST01000319.1 [span_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> 
- **GCA_000715385.1_Rhizoctonia_solani_123E_v1** Basidiomycota:PR AZST01000540.1 [core_model_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> ; models pheromone_receptor:438bp/44.44%->858bp/43.26%
- **GCA_000827255.1_Suilu1** Basidiomycota:HD KN835132.1 [locus_class_changed,span_changed,gene_set_changed,core_model_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; models HD1:456bp/36.24%->1950bp/67.08%;HD2:474bp/33.54%->1677bp/63.7%; genes HD1,HD2 -> HD1,HD2,MIP1,beta_fg
- **GCA_001542265.1_ASM154226v1** Basidiomycota:redPR LNKU01000018.1 [call_gained]: absent / -> called A2/high; margin  -> ; best score  -> 
- **GCA_001542265.1_ASM154226v1** Basidiomycota:redHD LNKU01000002.1 [call_gained]: absent / -> called undetermined/high; margin  -> ; best score  -> 
- **GCA_002092955.1_ASM209295v1** Basidiomycota:PR NCVV01000007.1 [span_changed,gene_set_changed,core_model_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> ; models fungal_mating_type_pheromone:added;pheromone_receptor:918bp/48.0%->1005bp/47.11%; genes caax_precursor,pheromone_receptor -> caax_precursor,fungal_mating_type_pheromone,pheromone_receptor
- **GCA_002900995.3_Gabo_G3** Basidiomycota:PR CM035310.1 [verification_changed,span_changed,gene_set_changed,core_model_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> ; models fungal_mating_type_pheromone:added;pheromone_receptor:951bp/52.4%->945bp/63.67%; genes caax_precursor,pheromone_receptor -> caax_precursor,fungal_mating_type_pheromone,pheromone_receptor
- **GCA_002900995.3_Gabo_G3** Basidiomycota:HD CM035305.1 [call_gained,span_changed,gene_set_changed]: withheld:modelled_gene_bar undetermined/ -> called undetermined/high; margin  -> ; best score  -> ; models no_gene_evidence_in_baseline; genes HD1,HD2 -> HD1,HD2,MIP1,beta_fg
- **GCA_002900995.3_Gabo_G3** Basidiomycota:HD CM035305.1 [call_gained,span_changed,gene_set_changed]: withheld:modelled_gene_bar undetermined/ -> called undetermined/high; margin  -> ; best score  -> ; models no_gene_evidence_in_baseline; genes HD1,HD2 -> HD1,HD2,MIP1,beta_fg
- **GCA_002995415.1_AZ12-009_Rhivin_1.0** Basidiomycota:HD PVYU01002681.1 [locus_class_changed,span_changed,gene_set_changed,core_model_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; models HD1:915bp/31.79%->1935bp/79.59%;HD2:462bp/37.01%->1695bp/71.99%; genes HD1,HD2 -> HD1,HD2,MIP1,beta_fg
- **GCA_004918335.1_ASM491833v1** Basidiomycota:wallMAT SPNZ01000038.1 [call_gained]: absent / -> called v1/high; margin  -> ; best score  -> 
- **GCA_019143615.1_ASM1914361v1** Basidiomycota:Aalpha JAGVSI010001513.1 [gene_set_changed,core_model_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; models HD1:99bp/51.52%->99bp/57.58%; genes HD1,HD2,Y,Z -> HD1,HD2,MIP1,Y,Z
- **GCA_019143615.1_ASM1914361v1** Basidiomycota:HD JAGVSI010001121.1 [span_changed,core_model_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; models HD1:462bp/33.55%->342bp/38.6%
- **GCA_019143615.1_ASM1914361v1** Basidiomycota:Balpha JAGVSI010000976.1 [gene_set_changed,core_model_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; models fungal_mating_type_pheromone:added; genes bap3-1,bap3-3,bar3,bbp2-1,bbp2-2,bbp2-5,bbp2-6,bbp2-8,bbp2_b,bbr2,caax_precursor,pheromone_receptor -> bap3-1,bap3-3,bar3,bbp2-1,bbp2-2,bbp2-5,bbp2-6,bbp2-8,bbp2_b,bbr2,caax_precursor,fungal_mating_type_pheromone,pheromone_receptor
- **GCA_022264835.1_Melame1** Basidiomycota:rustHD JAIFHY010000001.1 [call_gained]: absent / -> called undetermined/high; margin  -> ; best score  -> 
- **GCA_050403325.1_S_johnsonii_JCM1840v1** Basidiomycota:redHD BAABVI010000025.1 [call_gained]: absent / -> called undetermined/high; margin  -> ; best score  -> 
- **GCA_056247695.1_ASM5624769v1** Basidiomycota:HD JBVYWQ010000019.1 [locus_class_changed,span_changed,gene_set_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; genes HD1,HD2 -> HD1,HD2,MIP1,beta_fg
- **GCA_056247695.1_ASM5624769v1** Basidiomycota:PR JBVYWQ010001089.1 [call_gained,gene_set_changed]: withheld:below_fraction_floor undetermined/ -> called undetermined/medium; margin  -> ; best score  -> ; models no_gene_evidence_in_baseline; genes caax_precursor,pheromone_receptor -> caax_precursor,fungal_mating_type_pheromone,pheromone_receptor
- **GCA_943193715.1_gfAgaBisp1.1** Basidiomycota:HD OW971906.1 [locus_class_changed,span_changed,gene_set_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; genes HD1,HD2 -> HD1,HD2,MIP1,beta_fg
- **GCA_963514055.1_Cersu_MES15020_3.0_ref** Basidiomycota:PR CAUPSP010000007.1 [span_changed,gene_set_changed,core_model_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> ; models fungal_mating_type_pheromone:added;pheromone_receptor:966bp/56.92%->1062bp/67.24%; genes caax_precursor,pheromone_receptor -> caax_precursor,fungal_mating_type_pheromone,pheromone_receptor
- **GCA_963514055.1_Cersu_MES15020_3.0_ref** Basidiomycota:HD CAUPSP010000001.1 [call_gained,span_changed,gene_set_changed]: withheld:modelled_gene_bar undetermined/ -> called undetermined/high; margin  -> ; best score  -> ; models no_gene_evidence_in_baseline; genes HD1,HD2 -> HD1,HD2,MIP1,beta_fg
- **GCF_000149925.1_ASM14992v1** Basidiomycota:rustHD NW_003526561.1 [call_gained]: absent / -> called undetermined/high; margin  -> ; best score  -> 
- **GCF_000182895.1_CC3** Basidiomycota:Aalpha NW_003307544.1 [span_changed,gene_set_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; genes HD1,HD2,Y,Z -> HD1,HD2,MIP1,Y,Z,beta_fg
- **GCF_000218685.1_v1.0** Basidiomycota:PR NW_006763292.1 [span_changed,gene_set_changed,core_model_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> ; models fungal_mating_type_pheromone:added;pheromone_receptor:966bp/61.01%->1158bp/60.26%; genes caax_precursor,pheromone_receptor -> caax_precursor,fungal_mating_type_pheromone,pheromone_receptor
- **GCF_000218685.1_v1.0** Basidiomycota:PR NW_006763293.1 [core_model_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> ; models pheromone_receptor:579bp/47.4%->957bp/54.92%
- **GCF_000218685.1_v1.0** Basidiomycota:PR NW_006763304.1 [call_lost,span_changed]: called undetermined/medium -> withheld:below_fraction_floor undetermined/; margin  -> ; best score  -> ; models no_gene_evidence_in_candidate
- **GCF_000218685.1_v1.0** Basidiomycota:HD NW_006763290.1 [call_gained,span_changed,gene_set_changed]: withheld:modelled_gene_bar undetermined/ -> called undetermined/high; margin  -> ; best score  -> ; models no_gene_evidence_in_baseline; genes HD1,HD2 -> HD1,HD2,MIP1,beta_fg
- **GCF_000264905.1_Stehi1** Basidiomycota:PR NW_006763141.1 [call_lost,span_changed]: called undetermined/medium -> withheld:below_fraction_floor undetermined/; margin  -> ; best score  -> ; models no_gene_evidence_in_candidate
- **GCF_000264905.1_Stehi1** Basidiomycota:HD NW_006763131.1 [call_gained,span_changed,gene_set_changed]: withheld:modelled_gene_bar undetermined/ -> called undetermined/high; margin  -> ; best score  -> ; models no_gene_evidence_in_baseline; genes HD1,HD2 -> HD1,HD2,MIP1,beta_fg
- **GCF_000271585.1_Trametes_versicolor_v1.0** Basidiomycota:PR NW_007360328.1 [verification_changed,span_changed,gene_set_changed,core_model_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> ; models fungal_mating_type_pheromone:added;pheromone_receptor:954bp/54.46%->1452bp/100.0%; genes caax_precursor,pheromone_B43,pheromone_receptor -> caax_precursor,fungal_mating_type_pheromone,pheromone_B43,pheromone_receptor
- **GCF_000271585.1_Trametes_versicolor_v1.0** Basidiomycota:HD NW_007360321.1 [call_gained,span_changed,gene_set_changed]: withheld:modelled_gene_bar undetermined/ -> called undetermined/high; margin  -> ; best score  -> ; models no_gene_evidence_in_baseline; genes HD1,HD2 -> HD1,HD2,MIP1,beta_fg
- **GCF_000271605.1_Fomme1** Basidiomycota:HD NW_006760390.1 [locus_class_changed,span_changed,gene_set_changed,core_model_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; models HD1:252bp/32.14%->378bp/34.13%; genes HD1,HD2 -> HD1,HD2,MIP1
- **GCF_000292625.1_Dacryopinax_sp._DJM_731_SSP1_v1.0** Basidiomycota:PR NW_024467210.1 [span_changed,core_model_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> ; models pheromone_receptor:894bp/36.05%->912bp/37.46%
- **GCF_000320585.1_Heterobasidion_irregulare_v2.0** Basidiomycota:HD NW_009258197.1 [locus_class_changed,span_changed,gene_set_changed,core_model_changed]: called undetermined/high -> called undetermined/high; margin  -> ; best score  -> ; models HD1:600bp/33.5%->1950bp/87.46%;HD2:501bp/25.61%->1707bp/81.72%; genes HD1,HD2 -> HD1,HD2,MIP1,beta_fg
- **GCF_000320585.1_Heterobasidion_irregulare_v2.0** Basidiomycota:PR NW_009258203.1 [verification_changed,span_changed,gene_set_changed,core_model_changed]: called undetermined/medium -> called undetermined/medium; margin  -> ; best score  -> ; models fungal_mating_type_pheromone:added;pheromone_receptor:966bp/58.18%->1731bp/100.0%; genes caax_precursor,pheromone_receptor -> caax_precursor,fungal_mating_type_pheromone,pheromone_receptor
