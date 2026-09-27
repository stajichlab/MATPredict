# btbA in both idiomorphs, and no allele ignored for homothallic calls

2026-09-27. Code 4457c3a (PR #9): dab2eea (homothallic guard), 4457c3a (btbA).
Tests 790 pass.

## Tier replay on the real code (replay_real.py, replay/)

Reproduction unchanged (100% in 7 runs, Mucoromycota 241/245). Risers vs
1e69622: Serinales 39 -> 24; the 15 removed are all homothallic_candidate
(incl. GCA_030462985.1, now medium). 0 of 139 Serinales homothallic calls rise.
Other runs unchanged (cap6 10, Dothideomycetes 13 + 7, Wallemia 17).

## Mucoromycota scan, 4457c3a vs 1a00b0a (compare_vs_1a00b0a.tsv)

245 calls both; 0 new, 0 lost; 207 unchanged.
- 34 medium -> high (2 also partial_locus -> mat_locus): 33 Minus calls with a
  weak sexP cross-hit ignored (mostly Rhizopus with btbA), 1 Plus call
  (GCA_060309305.1) with sexM ignored. With btbA Plus-only, a sexM + btbA
  locus named both idiomorphs, so the rule could not apply.
- 4 idiomorph -> undetermined: all from the 25-bit classifier floor
  (533e266), not from btbA: Syncephalastrum racemosum NRRL 2496 (12.0),
  Phascolomyces articulosus (18.1), Benjaminiella poitrasii (one of three
  loci, 11.7), Mucor ardhlaengiktus CBS 210.80 (one locus, 24.2).

Zygo 23: scaffold 23/23 locus, 23/23 idiomorph; contig 23/23, 23/23.
