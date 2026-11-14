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
33 genomes searched (one of 34 was suppressed): Taphrina 13, Schizosaccharomyces 8 searched, 9 listed (not S. pombe, which has curated records),
Protomyces 8, Saitoella 2, Neolecta 1, Savitreella 1. Pneumocystis is not in `samples.csv` by class and was not run.

| result | genomes |
|---|---|
| `mat1` family (Schizosaccharomyces curated), mat1-Mc/Pc/Pi found | 7 of the 8 Schizosaccharomyces genomes searched; none in *S. lindneri* |
| `PM` family, matMc + matPi (curated from *Taphrina deformans*), idiomorph M, gene only | 4 Protomyces genomes (inouyei, macrosporus, two lactucae-debilis), *Taphrina confusa* |
| `MAT` family (Pezizomycotina references) | 16 calls across Taphrina, Protomyces, Saitoella, Neolecta |
| no call at all | 10 genomes: *Taphrina betulina*, *T. communis*, *T.* sp. SM11, *T. wiesneri* (one genome), *Protomyces arabidopsidicola*, *P. gravidus*, *P. pachydermus*, *Saitoella complicata*, *Savitreella*, *S. lindneri* |

Strength of the evidence (identity / coverage to the reference, from `detection_report.yaml`):
- `PM` calls: matMc 44 to 69% identity at 74% coverage; matPi 42 to 59% at 33 to 89% coverage. Strongest: *Taphrina confusa*
  (matMc 69.2% / 74.5%, matPi 42.2% / 89.0%). Confidence `medium`. Retracted after the structural follow-up below: the Protomyces calls are not distinguishable from background HMG genes; only the *T. confusa* matMc gene stays a lead.
- `MAT` family calls: identity 17 to 55%, coverage mostly 10 to 79%, references from Pezizomycotina. I do not treat them as
  mating-type calls. The 5 *T. deformans* genomes repeat one MAT1-2-1 hit at 44.9% identity and 17.9% coverage.
- `mat1` calls agree with the curated *S. pombe* family; they are not an independent check.

Not verified: that the `PM` hits contain the HMG domain, that the two genes are adjacent in the same layout as *T. deformans*,
and that the 10 empty genomes lack a locus (the detector recall for distant lineages is limited: see the holdout notes).

## Next
- (Done, see the follow-up section) Pfam and gene-order check of the `PM` leads (Protomyces x4, *T. confusa*), then a second look at the 10 empty genomes with the
  flank-first scan using matMc/matPi as controls.
- Draft candidate records only after a curator reads the structural check.

## Lineage coverage survey (added 2026-11-14)
BFD `samples.csv` (23,679 rows) was matched to the genomes in the flank-first scan (496). Scanned: Basidiobolomycota 18, Blastocladiomycota 14,
Chytridiomycota 83, Entomophthoromycota 9, Glomeromycota 59, Kickxellomycota 177, Mortierellomycota 100, Neocallimastigomycota 25, Zoopagomycota 11.
Not scanned before and now scanned with the same stage 1 (job 32090926, 141 genomes): Aphelidiomycota 2, Calcarisporiellomycota 1,
Monoblepharomycota 3, Olpidiomycota 1, Rozellomycota 1, one unclassified, Microsporidia 132. Not scanned because the BFD sample name did not map
to an annotation folder (name matching is imperfect, so counts are approximate): some Chytridiomycota, Entomophthoromycota, Kickxellomycota and
Neocallimastigomycota genomes, Sanchytriomycota 2, 31 Microsporidia and 2 Aphelidiomycota. Mucoromycota (Mucoromycotina) has curated MAT records and was not part of this survey.

Result table: `2026-11-14_mat-discovery-tables/gap_lineages_summary.tsv` (same columns as the lineage note of PR #351).
- Small phyla (9 genomes): no HD1+HD2 locus and no 2xHomeo locus. One `HMG+Homeo` locus each in the Calcarisporiella genome and the Rozella genome, and
  `sexM/sexP` control hits in 4 genomes. A single genome per lineage; not tested against chance.
- Microsporidia (132 genomes): 65 genomes have two Pfam Homeodomain genes within 5 genes. Random placement of the same number of genes gives a mean of
  18.8 genomes (maximum 28 in 200 shuffles), so the adjacency is more than chance. The median genome has 8 Homeodomain-hit genes, 79 of 154
  neighbouring pairs are on the same strand, and 104 pairs have no gene between them. That fits a tandem-duplicated family, not a MAT locus. No further test done.

## Allomyces probe (added 2026-11-14)
`scripts/flank_idiomorph_probe.py` (tests in `tests/test_flank_idiomorph_probe.py`). Stage 1 had found 35, 45 and 63 HMG/Homeodomain candidate genes in
*A. arbusculus* Burma 1F, *A. javanicus* California 12 and *A. macrogynus* ATCC 38327, all as single-gene loci (no adjacent pair). The new test asks, for each
candidate, whether its flanking genes form a syntenic block in the other genomes and whether the candidate is conserved there, replaced by another
HMG/Homeodomain gene, or missing.
- 84 candidate genes; 152 gene-by-genome comparisons were testable. Not conserved (absent or replaced): 56 of 152 = 0.37.
- Baseline, random non-candidate genes through the same test: 944 of 2,529 = 0.37. The candidates are not more often non-conserved than ordinary genes
  (per pair 0.19 to 0.48 against 0.13 to 0.54).
- 15 genes are absent or replaced in both other genomes, for example F5BE82F7_007902 (HMG_box, *A. macrogynus*; absent in both others, 4 flank genes each in a syntenic block).
  With a 0.37 background rate these are expected by chance; none is a call. F5BE82F7_000331 (*A. macrogynus*, HMG_box with a weak sexP hit, E 9.3e-8) was not
  singled out by the test.
- Not tested: gene order beyond 25 genes, divergent alleles that fall below diamond E 1e-5, and assembly gaps. No new Allomyces candidate. Allomyces results are negative.

## Structural follow-up of the PM leads and rescan of the 10 empty genomes (added 2026-11-14)
Correction: the `PM` calls are in 4 Protomyces genomes plus *T. confusa* (5 genomes), not 5 Protomyces genomes plus *T. confusa* as written above before this edit.

Script `scripts/pm_lead_followup.py` (tests in `tests/test_pm_lead_followup.py`); job 32091035; tables in `2026-11-14_mat-discovery-tables/`
(`pm_lead_summary.tsv`, `pm_lead_genes.tsv`, `empty_scan_genes.tsv`, `empty_scan_pairs.tsv`). Pfam is `hmg_homeo.hmm` (HMG_box, HMG_box_2, Homeodomain, GA cutoff).

PM leads (gene under each detector hit, tested with Pfam):
| genome | matMc gene | HMG-box | matPi gene | homeodomain |
|---|---|---|---|---|
| *Protomyces inouyei* | F7C0A940_005637 | yes | no gene under the hit | n/a |
| *P. lactucae-debilis* NRRL YB-4353 | F21A9957_003986 | yes | no gene under the hit | n/a |
| *P. lactucae-debilis* 12-1054 | F003EF48_006292 | yes | no gene under the hit | n/a |
| *P. macrosporus* | no gene under the hit | n/a | FCFD2032_003005 | no |
| *Taphrina confusa* | F57803F5_005977 | yes | no gene under the hit | n/a |
- In no genome do a matMc gene and a matPi gene sit next to each other. The two-gene `PM` structure is not confirmed anywhere. The matMc genes are HMG-box genes; the matPi hit has no Pfam homeodomain in *P. macrosporus* and no predicted gene in the other four.
- The Protomyces hits are not specific. The matMc hits are 44% identical, and the ten genomes with no call also have HMG genes that match the curated mat1-Mc at 46 to 49% (for example *P. gravidus* 46.1%, *P. pachydermus* 48.6%). The Protomyces `PM` calls are indistinguishable from background HMG genes in the same genomes.
- *T. confusa* is different: its matMc gene is 69.2% identical to the curated *T. deformans* matMc.

Rescan of the 10 genomes with no call (Pfam plus diamond, E 1e-3, against the 14 curated PM and mat1 proteins):
- Each genome has 5 to 7 HMG-box genes and 2 to 7 homeodomain genes (domain genes per genome 7 to 14). No genome has an HMG and a homeodomain gene, or two domain genes, within 5 genes of each other (`scan_pairs` is empty).
- *Taphrina wiesneri* CBS 275.28 has a strong matMc-like gene (72.5% identity, E 7.6e-110) that the detector did not report as a locus, because no second gene was found next to it. No other genome has a hit above 50% identity to matMc.
- Best matPi hits (31 to 55% identity) are generic homeodomain genes that every genome has; they carry no locus information.

Conclusion, with limits: no Taphrinomycotina genome tested shows a two-gene PM locus. A matMc-like gene (about 70% identity to *T. deformans*) is present in *T. confusa* and *T. wiesneri* CBS 275.28 with no partner nearby; the 3 other *Taphrina* genomes with no call have no such gene. I did not look for a matPi that sits far from matMc (the curated record places matMc and matPi on separate contigs of the old assembly), and the layout
test cannot see a locus split over contigs. The probe finds no new candidate to curate; the only lead worth keeping is the matMc-like gene in the two *Taphrina* genomes.
