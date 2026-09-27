# 007. Serinales-wide scan: homothallic_candidate calls

- **Category:** homothallism candidate
- **Status:** candidate
- **Lineage:** Ascomycota, Saccharomycotina, Serinales

## Summary
The Serinales-wide scan labels 139 calls `homothallic_candidate`. They are
excluded from the allele-absent confidence rule, so none is promoted by
ignoring one of its alleles. The per-species breakdown has not been written up.

## Evidence
- Count 139: `results/2026-09-27_btbA_homothallic/NOTE.md` ("0 of 139
  Serinales homothallic calls rise"), on the scan
  `results/2026-09-26_serinales_all_882aa01/` (main checkout).
- Species breakdown: unverified (not tabulated in any note).
- Debaryomyces caution: 42 `homothallic` genotypes had become
  `a+homothallic` through a spurious call at the PAP1-OBP1-PIK1 block, which
  is not at the MTL (entry 019).

## Method that found it
Rule `338ee18` (entry 005); guard `dab2eea` (homothallic calls never have an
allele ignored).

## Verification done / still open
Open: tabulate by species; separate true two-idiomorph loci from merged or
collapsed assemblies (entry 011 shows C. albicans assemblies collapse
heterozygous MTL, which argues for care in both directions).

## Limits
Gene content only; no per-species review yet.

## Related
Entries 004, 005, 011, 019.
