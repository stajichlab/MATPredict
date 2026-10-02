# Regression check: LCG Mucoromycotina 94d9d37 vs f9ad092 -- lcg

- Genomes compared: 621; loci rows: 3576; loci with a change: 6 (3 touch a call, 3 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 0 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| call_lost | 3 |

## Changes to calls

- **Fennellomyces_verticellatus_RSA_2442_Plus** Mucoromycota:MAT scaffold_26 [call_lost]: called undetermined/medium -> withheld:mat_gene_gate undetermined/; margin 15.4 -> 15.4; best score 44.2 -> 44.2; models no_gene_evidence_in_candidate
- **Fennellomyces_verticellatus_RSA_2446** Mucoromycota:MAT scaffold_287 [call_lost]: called undetermined/medium -> withheld:mat_gene_gate undetermined/; margin 15.4 -> 15.4; best score 44.2 -> 44.2; models no_gene_evidence_in_candidate
- **Thamnostylum_nigricans_RSA_1405-** Mucoromycota:MAT scaffold_274 [call_lost]: called Plus/medium -> withheld:mat_gene_gate Plus/; margin 48.6 -> 48.6; best score 76.3 -> 76.3; models no_gene_evidence_in_candidate
