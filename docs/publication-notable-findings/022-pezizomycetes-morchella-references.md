# 022. Morchella references and SLA2/APN2 flanks raise Pezizomycetes calls from 27 to 98

- **Category:** curation-first
- **Status:** verified (synteny)
- **Lineage:** Ascomycota, Pezizomycetes

## Summary
Two Morchella importuna MAT records and SLA2/APN2 as optional flanks raised
Pezizomycetes calls from 27 to 98 of 172 genomes, with no call lost. 89 of the
98 calls sit next to each genome's SLA2/APN2 ortholog, located independently.

## Evidence
- Records `1174673_ypl6-1_MATtub_MAT1-2`, `1174673_ypl6-3_MATtub_MAT1-1` (Chai
  et al. 2017, doi:10.1007/s11557-017-1309-x; deposits KY782629.1 / KY782630.1;
  7/7 genes at 100%).
- Discinaceae 1 -> 41 of 72; Morchellaceae 14 -> 43 of 60; Tuberaceae 12/12.
- Source: `docs/notes/2026-09-25_morchella-sla2-and-homothallic-candidates.md`;
  `results/2026-09-25_pezizomycetes_d689647/compare_vs_2026-09-24.txt`;
  synteny `results/2026-09-24_bar_synteny/`.

## Method that found it
Frozen worktree `run-d689647` (commits `81ea22a`, `ca619b6`).

## Verification done / still open
Done: independent SLA2/APN2 adjacency.

## Limits
Most new calls are medium; runtime 22 -> 60 s median per genome.

## Related
Entry 005 (Hydnotrya candidates from the same runs).
