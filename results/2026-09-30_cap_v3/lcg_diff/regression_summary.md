# Regression check: V3 cf5cb1e vs ad3ea9b -- LCG Mucoromycotina -- lcg

- Genomes compared: 621; loci rows: 3368; loci with a change: 13 (1 touch a call, 12 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 0 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| call_gained | 1 |
| classifier_input_changed | 1 |
| classifier_shift | 1 |

## Changes to calls

- **Mucor_pusillus_NRRL_A-13674** Mucoromycota:MAT scaffold_161 [call_gained,classifier_input_changed,classifier_shift]: withheld:modelled_gene_bar+polish_capped Plus/ -> called Plus/medium; margin 97.5 -> 90.5; best score 224.8 -> 200.1; models no_gene_evidence_in_baseline
