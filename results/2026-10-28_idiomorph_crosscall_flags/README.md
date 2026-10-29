# Genomes with loci from two Saccharomycotina MAT systems (2026-10-28)

`crosscall_genomes.tsv` (181 genomes; `scripts/idiomorph_crosscall_flags.py` on the committed v0.7.0 tables): genomes with loci called by two or more of the families MATsc, MTL, MATog, MATkp.
- **All 181 are `phylum_fallback` genomes**: no lineage record matched, so every Ascomycota family was tried. None of the lineage- or direct-routed genomes is in the list.
- **101 have only duplicate calls** (the secondary locus has the same contig, start and end as a locus of the main family): one locus called by two panels. `suggested_action` `duplicate_family_call` (78) or `duplicate_family_call+reroute` (23).
- **80 have a secondary locus at another position** (`curate`, 74; `curate+reroute`, 6): a second copy, a silent cassette or a second cross-call; needs a look.
- `+reroute`: a genus with its own family (Ogataea -> MATog, Komagataella -> MATkp) whose genome was routed by phylum fallback.

## Worked example: *Ogataea glucozyma* GCA_003706155.2 (routing `phylum_fallback`)
Both MATsc and MTL call contig `PPKO02000069.1:20-23514`, as MATalpha and as alpha, each `high` confidence, `idiomorph_gene_only`, with identical boundaries: MATALPHA1 / MATALPHA2 at 30.7% / 47.4% identity (MATsc references) and MTLalpha1 / MTLalpha2 at 32.9% / 45.3% (MTL references). The same locus is called twice, at 30-47% identity, by two panels that are not the right ones for *Ogataea* (the MATog family was not tried for this genome). A second MATsc locus (MATa, `partial_locus`, medium) lies on `PPKO02000123.1` next to `sla2`; the MTL panel's weak MTLA1 hit at the same contig is withheld. 40 more single-gene hits are withheld.

## Why it matters
- Counts by family and idiomorph, the website and any label double-count these loci.
- Pichiales (*Ogataea* 72 of 112 genomes, *Pichia* 157 of 158) are mostly phylum-fallback routed.

## Proposed
1. **Re-run** the 181 genomes with routing that tries MATog for *Ogataea* and the right family by genus for *Pichia* and relatives (or add a Pichiales-level route), and compare.
2. **Dedupe rule** in the report/table step: when two families call identical coordinates, keep the better-scoring family and mark the other `duplicate_family_call`.
3. **Curate** the 80 genomes with secondary loci elsewhere, starting with one *Ogataea* and one *Pichia* genome, and with a reference for *Pichia* (none exists).

## Dedupe rule (implemented, output only)
`scripts/dedupe_family_calls.py` marks loci that two families called at identical coordinates (`duplicate_family_call`, `kept_family`, `kept_idiomorph`) and keeps the better-scoring family: locus class, then confidence tier, polished genes, genes found, the family with more loci in the genome, then the family name. `loci_v070_ascomycota_deduped.tsv` is the v0.7.0 Ascomycota loci table with the three columns added (no row removed): 185 loci in 164 genomes are duplicates, by dropped family MTL 86, mat2 69 (S. pombe, kept as mat1), MATsc 27, MAT 2, MATtub 1. It is not wired into the site build or the clean-run tables; counts with and without it can be taken from the columns. Tests: `tests/test_dedupe_family_calls_script.py`.
