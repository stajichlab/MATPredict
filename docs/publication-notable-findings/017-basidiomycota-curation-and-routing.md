# 017. Basidiomycota: curated orders are called at 92-100%; first rust and red-yeast references

- **Category:** curation-first
- **Status:** verified (held-out genomes)
- **Lineage:** Basidiomycota

## Summary
In a 3,269-genome Basidiomycota run, orders routed to their own curated family
are called at 92-100%, while phylum-fallback orders mostly reach 0-20%. The
first Pucciniales and Sporidiobolales references lift those orders on genomes
that supplied no record.

## Evidence
- `results/2026-09-26_basidiomycota_full/ANALYSIS.md` (committed on
  basidio-anchors `efa9afa`): Agaricomycotina 83.0% of 2,417; Ustilaginomycotina
  95.6% of 315; Pucciniomycotina 13.2% of 486; Wallemiomycotina 0% of 51 (before
  entry 001). Curated orders: Agaricales 96.7%, Boletales 94.5%, Tremellales
  94.4%, Polyporales 97.5%, Russulales 93.5%, Ustilaginales 100%. 875 of 896
  uncalled fail the bar with 0 modelled genes (reference gap).
- `docs/notes/2026-09-26_pucciniomycotina-curation.md` (curation-puccinio):
  Sporidiobolales 4 -> 15 of 15; Pucciniales 0 -> 10 of 11, and 0 -> 7 of 8
  excluding record genomes. Records `redPR` (7, tier 1, Coelho et al. 2010,
  PMID 20700437; 2011, PMID 21880139) and `rustHD` (4, tier 2, Cuomo et al.
  2017, PMID 27913634).
- New Russulales, Boletales, Polyporales HD references (Heterobasidion
  KF280353.1, Rhizopogon AB646132.2, Phanerochaete HQ188438.1) and MIP1/beta-fg
  anchors: `docs/notes/2026-09-26_basidiomycota-anchors.md` (basidio-anchors).
  Agaricomycotina pilot 8/18 -> 13/18.

## Method that found it
Order-rank taxonomic routing; detection with frozen worktrees
(`run-ad1f865`, `run-293640d`).

## Verification done / still open
Open: pheromone-receptor loci outside Agaricales (entry 003); Microbotryales,
Boletales and Hymenochaetales receptors.

## Limits
Only 3 Basidiomycota:PR calls in the full run; receptor loci are called under
other family names (Balpha/Bbeta 106, aLocus 185, Tremellales MAT 385).

## Related
Entries 001, 002, 003, 009.
