# Flank-carried rule
Status: decided

## Results
- Serinales origin: 131/2,647 calls flank-carried; 52 Debaryomyces artefacts.
- Ascomycota audit (`results/2026-09-26_flank_rule_ascomycota/NOTE.md`): 89
  changed calls; two design faults (all-hits test, 3 kb window).
- Fix b0898d9: strongest core hit, per-family window (3 kb MTL, 20 kb others).
- Bitscore floor evaluation (`results/2026-09-27_flank_bitscore_floor/NOTE.md`):
  39 bits keeps 7/8 real loci, 0/52 known noise, 21/89 audit calls.
- Implemented 51ec961 (`results/2026-09-27_flank_bitscore_implemented/NOTE.md`):
  21/89 kept as predicted; 12 newly kept Ascomycota calls all real on gene order.

## Limits
Floor tuned on 8 real loci, 5 of them the Umbelopsis loci it rescued.
