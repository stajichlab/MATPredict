# Bd cross-isolate depth scan for a mating-type idiomorph (2026-11-14)

Status: negative result. No candidate idiomorph was found. Nothing here is curated.

Script: `scripts/bd_idiomorph_scan.py` (tests in `tests/test_bd_coverage.py`). Input: the 10 kb window tables of 58 (CLFT044
reference) and 59 (JEL423 reference) isolates with median depth of at least 20x. Genes come from the BFD funannotate GFF3
(scaffold IDs identical to the references used). Tables: `tables/idiomorph_regions_<ref>.tsv`, `tables/idiomorph_null_<ref>.tsv`.

## Method
Each window ratio is divided by its sequence median, which removes whole-sequence copy-number change. A window is FULL
(0.8 to 1.25), HALF (0.35 to 0.65) or LOW (below 0.25) per isolate. A candidate window has at least 5 FULL isolates and at least
3 LOW or HALF isolates. Neighbouring candidate windows merge; regions of 20 kb or more are kept. Expected signature of a
diploid idiomorph: FULL, HALF and LOW isolates in the same window, with few het SNPs in HALF isolates (only one copy aligns).

## Result
- CLFT044 coordinates: 71 regions. JEL423 coordinates: 70 regions.
- Most regions are HALF in many isolates, with a median of 20 to 60 het SNPs per 10 kb window in those isolates. That fits
  divergent sequence that maps poorly onto a reference of another lineage, not a hemizygous region.
- Regions with n_half >= 3 and median het in HALF isolates <= 3: 3 against CLFT044 and 5 against JEL423. None has a gene whose
  product or note matches HMG, homeodomain, mating, sex or pheromone. HALF isolates in these regions come from both CLFT-type and
  JEL-type groups, so lineage does not explain them. I did not examine their gene content further.
- The rotation null (each isolate rotated by a random offset, 100 permutations) gave a mean of 120 (CLFT044) and 110 (JEL423)
  regions, more than observed. The region count is therefore not a usable test statistic, and I make no significance claim.

## Limits
- No evidence that Bd has a two-idiomorph system. A scan like this assumes a presence/absence or hemizygous idiomorph.
- 10 kb windows. A small locus can be missed. Smaller windows need a new alignment, because no BAMs were kept.
- Isolates below 20x were dropped. Two of the references are from different lineages, which drives the HALF windows.
