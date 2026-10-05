# Regression check: regression run-tiered vs run-af3c615 -- record_sources

- Genomes compared: 23; loci rows: 938; loci with a change: 3 (1 touch a call, 2 are withheld on both sides).
- Classifier margin shifts below 5.0 bits are not listed as changes; 0 unchanged loci carry such a shift (see `margin_delta` in the TSV).
- Every withheld-only change is listed in `regression_withheld_changes.md`.

## Change types (loci that touch a call)

| change | loci |
|---|---|
| call_gained | 1 |

## Changes to calls

- **GCF_000011425.1_ASM1142v1** Ascomycota:MAT NC_066262.1 [call_gained]: withheld:modelled_gene_bar+polish_capped MAT1-1/ -> called MAT1-1/medium; margin  -> ; best score  -> ; models no_gene_evidence_in_baseline
