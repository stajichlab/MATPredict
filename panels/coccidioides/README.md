# Coccidioides panels

- `cocci_MAT1-1.fasta`, `cocci_MAT1-2.fasta` (use this pair): the *C. immitis* H538.4 (EF472259.1, offset 5,300 to the end, 8,905 bp) and RS (NW_004504310.1 record, offset 3,895 to the end, 8,362 bp) idiomorph regions, trimmed to the idiomorph (COX13/APN2/MAT1-1-1 and SLA2-side flank to CIMG_00407/MAT1-2-1/APN2/COX13) plus the 1.7 kb SLA2-side flank that both records share (the shared-depth control). The untrimmed records carried several kb of flank that exists in both idiomorph backgrounds but only in one record, which made the other idiomorph's k-mers look present (the `trace_` flags in the pilot).
- `immitis_*` (full records, used in the pilot), `posadasii_*` (separate GenBank deposits with no shared flank; gives `low_depth`; not for use).
