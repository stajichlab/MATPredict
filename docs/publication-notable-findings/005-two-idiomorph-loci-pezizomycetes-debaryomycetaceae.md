# 005. Loci carrying both idiomorphs' genes in Hydnotrya and Debaryomycetaceae

- **Category:** homothallism candidate
- **Status:** candidate
- **Lineage:** Ascomycota: Pezizomycetes (Discinaceae, Hydnotrya); Saccharomycotina (Debaryomycetaceae)

## Summary
Genome-only runs found loci where the two opposite-idiomorph core genes are
unrelated proteins (alpha box vs HMG box, or homeodomain vs alpha box), so they
cannot be one region seen twice. After a stricter rule (full-length models, no
cross-match), six genomes are labelled homothallic candidates.

## Evidence
Source: `docs/notes/2026-09-25_morchella-sla2-and-homothallic-candidates.md`;
`results/2026-09-25_homothallic_check/`.

| genome | species | genes at one locus |
|---|---|---|
| GCA_040803765.1 | Hydnotrya cerebriformis | MAT1-2-1 (2,578-3,156) + MAT1-1-1 (5,139-6,075), exonerate models, ~49% |
| GCA_040803455.1 | Hydnotrya variiformis | MAT1-1-1 (8,931-9,867) + MAT1-2-1 (11,840-12,418), ~49% |
| GCA_030564605.1 | Debaryomyces coudertii NRRL Y-7425 | MTLA1 + MTLA2 + MTLalpha1 |
| GCA_056149625.1 | Debaryomyces hansenii Wch | MTLA1 + MTLA2 + MTLalpha1 |
| GCA_046867635.1 | Schwanniomyces etchellsii CBS 6823 | MTLA2 + MTLalpha1 |
| (CBS767) | Debaryomyces hansenii | entry 004 |

Below the coverage floor (not labelled): Priceomyces fermenticarens
GCA_030574455.1, P. melissophilus GCA_030574515.1.

## Method that found it
Rule `338ee18`: two core genes of different gene_class, both full-length
models (>= 50% of their reference), no cross-match. A domain-only relaxation
would have labelled 113 loci in heterothallic-rich panels; the full rule
labels 0 of those (10 re-run genomes confirmed plain `mat_locus`).

## Verification done / still open
- Done: the false-positive check above.
- Open: mating tests or literature per species. Dirks et al. 2025 report one
  colocalised MAT1-1/MAT1-2 case in Discinaceae needing confirmation (per the
  2026-09-25 note).

## Limits
Gene content only. Whether these genomes self-mate is not shown.

## Related
Entries 004, 006, 007.
