# 004. Debaryomyces hansenii CBS767 carries a1, a2 and alpha1 at one MTL locus

- **Category:** homothallism candidate; assembly-or-annotation artefact (the record error that hid it)
- **Status:** candidate
- **Lineage:** Ascomycota, Saccharomycotina, Serinales (Debaryomycetaceae)

## Summary
The curated record for D. hansenii CBS767 had stated that no MTLalpha1 gene is
present. XP_460134.1 sits about 1.5 kb from MTLa1, inside the MTL locus, and is
MTLalpha1. The strain therefore carries genes of both idiomorphs at one locus.

## Evidence
- XP_460134.1 at `NC_006047.2:1,587,954-1,588,592`, about 1.5 kb from MTLa1;
  the only PF04769 (alpha box) protein in the proteome; reciprocal best hit of
  C. albicans MTLalpha1 (E=2.8e-25); CTG-clade placement in a gene tree.
  Source: annotation report A1; `results/2026-09-24_dhansenii_alpha1/EVIDENCE.md`.
- Record `4959_cbs767_MTL_A` record_version 2 (commit `ca619b6`) adds
  MTLalpha1 and `system: homothallic` as a gene-content candidate.
- Detection labels it `homothallic_candidate` under rule `338ee18`
  (`docs/notes/2026-09-25_morchella-sla2-and-homothallic-candidates.md`).

## Method that found it
Reciprocal best hit, Pfam domain search and gene tree
(`results/2026-09-24_dhansenii_alpha1/`).

## Verification done / still open
- Done: RBH, domain, tree, position.
- Open: whether the strain self-mates. Krassowski et al. 2019 describe the
  contiguous a+alpha arrangement in Debaryomyces (per the 2026-09-25 note).

## Limits
Gene content is not proof of homothallism.

## Related
Entries 005, 007; annotation report B3.1 (XP_460134.1 named only "DEHA2E19096p").
