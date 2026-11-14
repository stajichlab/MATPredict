# MAT discovery: search strategy and the Taphrinomycotina probe (2026-11-14)

Status: exploratory. Nothing is curated or accepted. Source: detector run, job 32090545, code commit 9c4759269f (worktree
`run-9c47592`), preset v0.9.1, `--withheld-detail-policy on`. Results: `results/2026-11-14_taphrinomycotina_probe/` (HPCC, not committed).

## Goal (user, 2026-11-14)
Find mating-type loci in lineages with none defined. Chytrids first (Bd the main target, with the depth/heterozygosity comparison),
Allomyces, then other lineages not yet searched. Do not repeat taxa already found to have no good candidate. Probe Taphrinomycotina.

## What is already done (not repeated)
- Flank-first scan of 496 BFD genomes (PR #351): Mortierellomycota none; Blastocladiomycota (14 genomes, 3 Allomyces) nothing above
  chance; Zoopagomycota, Kickxellomycota, Basidiobolus and Chytridiomycota give weak or lineage-specific hits only.
  The Neocallimastigomycota (25) and Entomophthoromycota genome lists were part of the scan, but the lineage note does not report them.
- Bd isolate depth/heterozygosity, parent of origin and an idiomorph depth scan (PR #352): no idiomorph candidate.

## Search strategy (proposal; the detector family gate is not changed)
1. **Detector probe with curated families** on genomes of a lineage close to a curated one. Reads out family, identity, coverage.
   Treat low identity or low coverage hits as leads only.
2. **Verify a lead from structure, not score**: the HMG or homeodomain gene must have the domain (Pfam), the genes must be adjacent
   as in the curated record, and a second genome of the same genus should show the same layout.
3. **Population signal where isolates exist**: for a species with many isolates (Bd), look for presence/absence or hemizygous
   windows across isolates (done for Bd; negative).
4. **Lineages with two or three genomes (Allomyces)**: no population signal is possible. Use the flank-first scan plus a search with
   pheromone-receptor and HMG-box families, then test layout across the three genomes. I have not run a new Allomyces search.
5. **Record unsearched lineages before scanning**: compare the BFD phylum/class list with the lists already scanned, scan only the gaps.

## Taphrinomycotina probe
33 genomes searched (one of 34 was suppressed): Taphrina 13, Schizosaccharomyces 9 (not S. pombe, which has curated records),
Protomyces 8, Saitoella 2, Neolecta 1, Savitreella 1. Pneumocystis is not in `samples.csv` by class and was not run.

| result | genomes |
|---|---|
| `mat1` family (Schizosaccharomyces curated), mat1-Mc/Pc/Pi found | 8 of 9 Schizosaccharomyces; none in *S. lindneri* |
| `PM` family, matMc + matPi (curated from *Taphrina deformans*), idiomorph M, gene only | 5 Protomyces genomes, *Taphrina confusa* |
| `MAT` family (Pezizomycotina references) | 16 calls across Taphrina, Protomyces, Saitoella, Neolecta |
| no call at all | 14 genomes, including *Taphrina betulina*, *T. communis*, *T.* sp. SM11, *T. wiesneri* (one genome), *Protomyces arabidopsidicola*, *P. gravidus*, *P. pachydermus*, *Saitoella complicata*, *Savitreella*, *S. lindneri* |

Strength of the evidence (identity / coverage to the reference, from `detection_report.yaml`):
- `PM` calls: matMc 44 to 69% identity at 74% coverage; matPi 42 to 59% at 33 to 89% coverage. Strongest: *Taphrina confusa*
  (matMc 69.2% / 74.5%, matPi 42.2% / 89.0%). Confidence `medium`. These are the leads worth following.
- `MAT` family calls: identity 17 to 55%, coverage mostly 10 to 79%, references from Pezizomycotina. I do not treat them as
  mating-type calls. The 5 *T. deformans* genomes repeat one MAT1-2-1 hit at 44.9% identity and 17.9% coverage.
- `mat1` calls agree with the curated *S. pombe* family; they are not an independent check.

Not verified: that the `PM` hits contain the HMG domain, that the two genes are adjacent in the same layout as *T. deformans*,
and that the 14 empty genomes lack a locus (the detector recall for distant lineages is limited: see the holdout notes).

## Next
- Pfam and gene-order check of the `PM` leads (Protomyces x5, *T. confusa*), then a second look at the 14 empty genomes with the
  flank-first scan using matMc/matPi as controls.
- Draft candidate records only after a curator reads the structural check.
