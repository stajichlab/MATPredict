# 015. Mortierellomycota and Kickxellomycota lack the Mucorales MAT arrangement

- **Category:** biology
- **Status:** candidate
- **Lineage:** Mortierellomycota, Kickxellomycota (Mucoromycota as control)

## Summary
The Mucorales tptA-HMG-rnhA arrangement appears in 0 of 100 Mortierellomycota
and 0 of 190 Kickxellomycota genomes, although the flank orthologs are present.
None of their 13 calls is supported, and none falls in the sexM or sexP clade.

## Evidence
- `results/2026-09-26_flank_synteny_ED/NOTE.md` (committed `c6ebcda`): tptA+rnhA
  pair / pair+HMG: Mucoromycota 139 / 117 of 293; Mortierellomycota 0 / 0 of
  100 (both orthologs confirmed in all 100); Kickxellomycota 0 / 0 of 190;
  chytrids 0 / 0 of 19.
- Tree: none of the 13 calls in either clade (`results/2026-09-26_sexMP_phylogeny/NOTE.md`).
- Calls labelled `verification: unverified` (`f737863`).

## Method that found it
Reciprocal best hit of flank genes against three Mucorales proteomes;
co-location test at 100/300 kb.

## Verification done / still open
Open: a relocated MAT locus cannot be excluded by this test.

## Limits
The test also fails in some Mucorales families (Lichtheimiaceae 0/30,
Umbelopsidaceae 0/14, Syncephalastraceae 0/9).

## Related
Entry 014.
