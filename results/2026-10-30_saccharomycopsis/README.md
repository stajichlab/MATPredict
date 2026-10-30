# Saccharomycopsis schoenii: MATss family and two candidate records (2026-10-30)

Curator decisions (2026-10-30): draft records from our own gene models (tier 3); a scope-only family per genus, as for MATkp, MATog and MATpk.

- **Family** `MATss` (scope: genus *Saccharomycopsis*, taxid 4943; scope-only, never searched on a phylum fallback) in `db/Ascomycota/order.yml`.
- **Records** (candidates, `needs_review`): `44080_gdd12_MATss_MATalpha` (left sub-cluster: MATalpha1, MATalpha2, sla2) and `44080_gdd12_MATss_MATa` (right palindromic cassette: sla2, MATa1). Homothallic. Strain GDD12, assembly GCA_054953955.1, chromosome IV CM143129.1.
- **Evidence:** Kriti et al., G3 2026 ([10.1093/g3journal/jkag067](https://doi.org/10.1093/g3journal/jkag067), PMID 41849659) for the SLA2-MAT architecture and locus spans; Dost et al., Microbiol Res 2024 ([10.1016/j.micres.2024.127691](https://doi.org/10.1016/j.micres.2024.127691), PMID 38492364) for homothallism. **Gene coordinates are ours**: the authors' GTF is not public. Models, method and limits: `../2026-10-29_ss_records/README.md`.
- **Genetic code:** NCBI table 12 (CUG = Ser); every gene carries `transl_table: 12`.
- **Checks:** `make_records.py` re-extracts every CDS from the window fasta (stop codon, no internal stop, expected length); `curate-db propose` passed; `curate-db validate` under Slurm (`validate.slurm`, job 31445666) gave `taxonomy_current` true. No protein accessions, so no per-gene sequence match.

## Open points for the curator
- Gene naming: both homeodomain ORFs are closest to MATa1-type references (phmmer 50 and 47 bits to *C. albicans* MTLA1); the left one is named MATalpha2 after the paper by position only. An HMG-like ORF (470,349-470,951) is not assigned. No MATa2 in the right cassette.
- Both SLA2 ORFs lack a start codon (5-prime partial). The mirrored arm of the palindrome is not recorded.
- The label registry has no entry for `Ascomycota:MATss`; a system code needs your approval before the records are labelled (frozen codes).
- Scope is the genus; *S. lassenensis* is heterothallic (Dost et al.), so narrow to the species (44080) if a replay over-calls.

## The other two genera (same session)
`MATpf` (*Petasospora fabianii*, species scope 36022; heterothallic, MATa JOY008 and MATalpha TBRC 4498; AUGUSTUS and miniprot models, no RNA-seq exists) and `MATbz` (*Barnettozyma*, genus scope 599802; one combined PB4207 record, MATalpha2 not modelled) were annotated in `../2026-10-30_petasospora/` and `../2026-10-30_barnettozyma/`. All five candidate records pass `curate-db propose` and `validate` (Slurm job 31446171: `taxonomy_current` true, sequence match not applicable because no record has a protein accession). All are tier 3 (our own models), `needs_review`.
