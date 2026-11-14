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
- In no genome do a matMc gene and a matPi gene sit next to each other. The two-gene `PM` structure is not confirmed anywhere. The matMc genes are HMG-box genes; the matPi hit has no Pfam homeodomain in *P. macrosporus* and no predicted gene in the other four. **[RETRACTED 2026-11-14, see the genome-level matPi search below: matPi is unannotated, not absent; a homeodomain region lies 7 to 14 kb from the matMc-like gene in 4 of these genomes]**
- The Protomyces hits are not specific. The matMc hits are 44% identical, and the ten genomes with no call also have HMG genes that match the curated mat1-Mc at 46 to 49% (for example *P. gravidus* 46.1%, *P. pachydermus* 48.6%). The Protomyces `PM` calls are indistinguishable from background HMG genes in the same genomes. **[RETRACTED 2026-11-14, see the genome-level matPi search below: the HMG gene is generic, but its adjacency to a homeodomain region is not background]**
- *T. confusa* is different: its matMc gene is 69.2% identical to the curated *T. deformans* matMc.

Rescan of the 10 genomes with no call (Pfam plus diamond, E 1e-3, against the 14 curated PM and mat1 proteins):
- Each genome has 5 to 7 HMG-box genes and 2 to 7 homeodomain genes (domain genes per genome 7 to 14). No genome has an HMG and a homeodomain gene, or two domain genes, within 5 genes of each other (`scan_pairs` is empty).
- *Taphrina wiesneri* CBS 275.28 has a strong matMc-like gene (72.5% identity, E 7.6e-110) that the detector did not report as a locus, because no second gene was found next to it. No other genome has a hit above 50% identity to matMc.
- Best matPi hits (31 to 55% identity) are generic homeodomain genes that every genome has; they carry no locus information.

Conclusion, with limits: no Taphrinomycotina genome tested shows a two-gene PM locus. A matMc-like gene (about 70% identity to *T. deformans*) is present in *T. confusa* and *T. wiesneri* CBS 275.28 with no partner nearby; the 3 other *Taphrina* genomes with no call have no such gene. I did not look for a matPi that sits far from matMc (the curated record places matMc and matPi on separate contigs of the old assembly), and the layout **[RETRACTED 2026-11-14, see the genome-level matPi search below: a matMc-like HMG gene with a homeodomain region next to it exists in 6 genomes]**
test cannot see a locus split over contigs. The probe finds no new candidate to curate; the only lead worth keeping is the matMc-like gene in the two *Taphrina* genomes. **[RETRACTED 2026-11-14, see the genome-level matPi search below: the leads are the HMG + homeodomain pairs below]**

## matPi search and positive control (added 2026-11-14)
Scan of the 5 PM-lead genomes, the 10 no-call genomes and the 4 *T. deformans* genomes for homeodomain genes similar to the curated matPi (job 32091046, 32091068; tables `pm_genomes_scan_genes.tsv`, `taphrina_deformans_scan_genes.tsv`, `empty_scan_genes.tsv`).
- Positive control: in all 4 *T. deformans* genomes the curated matMc matches an HMG-box gene at 100% identity. The curated matPi matches no gene better than 49.2%, a generic homeodomain gene on a different contig from the matMc gene. The curated *T. deformans* matPi has no exact match in any *T. deformans* genome. The record already says matPi is unresolved; this is a second reason to review it. I did not trace where the curated matPi sequence comes from.
- Every *Taphrina* genome has one homeodomain gene at 47.6 to 50.0% identity to the curated matPi, including the 3 genomes with no matMc-like gene. *T. confusa* (50.0%) and *T. wiesneri* CBS 275.28 (50.0%) have it on a different contig from the matMc-like gene, as in *T. deformans*. Identity cannot tell a matPi ortholog from the other homeodomain genes (the top hit stays within 2 points of the *T. deformans* value), so I cannot say whether these genes are matPi. No homeodomain gene is closer than 55.2% to the curated matPi in any genome tested (*Protomyces arabidopsidicola* 55.2%).
- Protomyces genomes with a PM call have no homeodomain gene above the scan threshold for matPi (E 1e-3), apart from the *P. arabidopsidicola* hit above.

## Follow-up of the leads in the other lineages (added 2026-11-14)
Chance test (`scripts/adjacent_pair_null.py`, 500 random placements per genome; `gap_adjacent_pair_null.tsv`). With 6 lineages and 3 pair types, 18 tests were run: one p of 0.05 or below is expected by chance.

| lineage | pair type | observed genomes | expected mean (max) | p |
|---|---|---|---|---|
| Microsporidia (132) | 2xHomeo | 65 | 19.0 (30) | 0.002 (floor of 500 draws) |
| Microsporidia | HMG+Homeo | 0 | 4.5 (10) | 1.0 |
| Monoblepharomycota (3) | 2xHMG | 2 | 0.2 (2) | 0.014 |
| Rozellomycota clade (1) | HMG+Homeo | 1 | 0.0 (1) | 0.050 |
| Calcarisporiellomycota (1) | HMG+Homeo | 1 | 0.5 (1) | 0.46 |
| Aphelidiomycota (2), Olpidiomycota (1) | all | 0 | 0.0 to 0.0 | 1.0 |

- **Microsporidia 2xHomeo** is a real excess (65 against 19.0). 55.4% of the pairs are on the same strand and 104 of 154 neighbours have no gene between them. This fits a tandem-duplicated family and not one locus per genome. Not followed further.
- **Monoblepharomycota 2xHMG** comes from *Hyaloraphidium curvatum* (one adjacent same-strand HMG pair) and *Gonapodya prolifera* (4 separate HMG pairs on 4 contigs). Four pairs in one genome point to a repeated family, not one MAT locus. All genes are annotated as hypothetical proteins. The other *Gonapodya* genome has no pair.
- **Rozellomycota clade (*Paramicrosporidium saccamoebae* KSL3)**: HMG-box gene FE6F46FC_001951 and homeodomain gene FE6F46FC_001954 on MTSL01000117.1 with two genes between them (the HMG gene on `-`, the homeodomain gene on `+`, divergent). Neither has a hit to the AMF HD1/HD2 or the Mucorales sexM/sexP controls at E 1e-5. This is the only compact HMG+homeodomain arrangement among the new lineages; it is one genome, uncorrected p = 0.05, and cannot be tested across genomes. Not a candidate record.
- **Calcarisporiella**: pair not above chance.
- **Aphelidiomycota, Olpidiomycota**: no pair.

## Strategy status (checked against the five steps above)
| step | status |
|---|---|
| 1. Detector probe with curated families | done for Taphrinomycotina (33 genomes) and Antarctolithica earlier |
| 2. Structural verification of leads | done for PM leads (Pfam, gene order, positive control) and for the new-lineage leads (gene order, products, chance test). Result: no lead survives except the matMc-like gene in two *Taphrina* genomes and the single *Paramicrosporidium* pair |
| 3. Population signal | done for Bd (negative) |
| 4. Small-genome-count lineages | done for Allomyces (negative against a random-gene baseline) |
| 5. Record unsearched lineages | done to the limit set by the data (see the coverage recount at the end); 50 genomes remain unscanned and 4 more have gene sets too small to test |
Open: a curator check of the matPi reference of the *T. deformans* record; any test of the *Paramicrosporidium* pair needs a second genome of that clade.

## Genomes that were not in any scan, and what was done (added 2026-11-14)
The first survey left 65 genomes (target phyla only) outside the 496 + 141 scanned genomes. A name-matching fix and a look at each folder (table `unscanned_genomes_status.tsv`) gave (65 genomes: 3 annotated + 48 skipped + 10 without folder + 4 without proteome; 58 stay unscanned):

| status | genomes | action |
|---|---|---|
| annotated, scannable | 3 | stage 1 run (job 32091100): *Homolaphlyctis polyrhiza*, *Vavraia culicis*, *Rhizophagus irregularis* DAOM_181602. Only the Rhizophagus genome has an HD1+HD2 locus: the known AMF locus, so it works as a positive control |
| BFD annotation skipped (`skipped_too_small`) | 48 | 4 assemblies (three of 11 to 24 Mb in 44 to 495 contigs, and a 1.9 Mb *Synchytrium endobioticum* assembly in 203 contigs) got Augustus gene prediction and stage 1 (below). The other 44 are fragments; they include *Globomyces pollinis-pini* Arg68 (21.4 Mb in 983 contigs, about 22 kb per contig), which I did not predict of hundreds to thousands of contigs (13 Kickxellomycota at about 8 Mb in 1,600 to 2,000 contigs; Microsporidia 0.7 to 7.5 Mb in 278 to 2,245 contigs; Sanchytriomycota 4.2 and 6.3 Mb in 549 and 897 contigs) and were not scanned. A flank-first scan needs several genes per contig |
| no annotation folder | 10 (Chytridiomycota 4, Entomophthoromycota 2, Glomeromycota 2, Microsporidia 1, one unclassified) | not scanned: no gene models. Includes a *Batrachochytrium dendrobatidis* genome, *Entomophthora muscae*, *Massospora cicadina*, *Gigaspora margarita*, *Mitosporidium daphniae* |
| folder without proteome | 4 (Entomophthoromycota 1, Microsporidia 3) | not scanned |

Augustus gene prediction (`scripts/augustus_gene_predict.slurm`, `scripts/augustus_proteins.py`, model `rhizopus_oryzae`, not trained on these lineages; jobs 32092338, stage 1 job 32092954; table `augustus_genomes_candidates.tsv`):
| genome | contigs | genes predicted | HMG/Homeodomain candidates | result |
|---|---|---|---|---|
| *Amoeboaphelidium protococcarum* FD95 | 44 | 4,898 | 4 | no HD1+HD2, no HMG+flank, no pair |
| *A. protococcarum* X5 | 155 | 9,397 | 5 | same |
| *Rozella allomycis* CSF55 | 495 | 3,527 | 9 | same |
| *Synchytrium endobioticum* 18_T1 | 203 (1.9 Mb assembly) | 821 | 2 | partial assembly; too little sequence to test |
The summary count for the four genomes: HD1+HD2 0, HMG+flank 0, 2xHomeo 0, HMG+Homeo 0, sexM/sexP control hits in 2 genomes. The total candidate counts in stage 1 (domain hits plus control hits; 11, 17, 15 and 2) are lower than in BFD-annotated genomes of related lineages (13 to 71: *Olpidium bornovanus* 13, *Aphelidium insulamus* 19, *Gonapodya* 31 and 47, *Calcarisporiella thermophila* 71). That is what a mismatched Augustus model would give, but I did not test it. It suggests that these predicted gene sets miss domain genes, so a negative here is weak evidence. The *Synchytrium* assembly is only 1.9 Mb, so it says nothing about the species. A proper test needs BFD-style annotation of these assemblies (funannotate) or RNA-based gene models.

Coverage conclusion: of the genomes in the target phyla that I could map, every one with gene models has been through at least one stage-1 scan, apart from those listed as not scanned above (58). The lineage conclusions stand only for those. In particular, for Sanchytriomycota, Kickxellomycota fragments, Microsporidia fragments and the Entomophthoromycota without annotation, nothing was tested.

## Coverage recount and final status (added 2026-11-14; supersedes the counts of the two sections above)
Recount by locus tag (the `LOCUSTAG` column of `samples.csv` against the gene IDs in all stage-1 outputs), which does not depend on folder names. The earlier name-based counts (58 unscanned; 10 genomes without an annotation folder) were partly a naming artefact: for example *Entomophthora muscae*, *Massospora cicadina*, *Gigaspora margarita* and a *Chytridiomycetes* sp. had been scanned in the first run under another folder name.
- Target phyla in `samples.csv`: 699 samples. After the first two scans, 61 had no stage-1 output by locus tag; 4 of these were the Augustus genomes of the previous section, so 57 were truly unscanned. All 57 lack BFD gene models (`unscanned_by_locustag.tsv`). Every target-phylum genome that has BFD gene models has been through stage 1.
- Contiguity of the 57 (`unscanned_assembly_stats.tsv`): N50 of 50 kb or more in 5, 20 kb or more in 9; one has no genome file.
- Seven more assemblies were given Augustus gene models (same caveat: `rhizopus_oryzae` model, not trained on these lineages) and stage 1: both *Globomyces* assemblies (chytrids), two *Mitosporidium daphniae* assemblies and three *Spraguea lophii* assemblies (jobs 32093014, 32093248; stage 1 jobs 32093239, 32093330).

| genome | contigs / N50 | genes | candidates | result |
|---|---|---|---|---|
| *Globomyces pollinis-pini* Arg68 | 983 / 51 kb | 8,406 | 17 HMG, 3 Homeodomain | 2xHMG pairs, conserved HMG cluster (below) |
| *Globomyces* sp. JEL0801 | 800 / 46 kb | 5,960 | 10 HMG, 2 Homeodomain | 2xHMG pairs, same cluster |
| *Mitosporidium daphniae* (ASM5257605) | 10 / 864 kb | 2,420 | 4 Homeodomain | no HMG gene, no pair |
| *M. daphniae* UGP1 | 273 / 35 kb | 2,359 | 7 | no pair |
| *Spraguea lophii* x3 | 143 to 158 / 97 to 99 kb | 356 to 362 | 0 | gene sets far too small (a microsporidian of 5.6 Mb should have about 2,000 genes); not testable |

Globomyces (chytrid) follow-up, with the one genus where two genomes can be compared (`flank_idiomorph_probe.py --ann`, tables `globomyces_probe.tsv`, `globomyces_baseline.tsv`):
- 32 HMG/Homeodomain candidates; 21 testable; 6 not conserved (0.29) against 112 of 728 random genes (0.15). Six against about 3 expected; I did not compute an exact p (a rough binomial tail is about 0.1), so this is not significant.
- The 2xHMG excess (2 of 3 genomes, p 0.008 in `adjacent_pair_null_aug2.tsv`) rests on 3 genomes and several tests and is not taken as evidence.
- The conserved block with an HMG cluster (Glopol1 MU621598.1 genes 7540 to 7545; JEL0801 JADGIX010000364.1 genes 5331 to 5335) holds two HMG gene types. The long type (263 aa) is 99.2% identical between the species and in the same position. The short type (138 aa) is 76.8% identical between the species and occurs as two copies 97.8% identical to each other in Glopol1 and one copy in JEL0801: a species-specific tandem duplication. Both species carry both genes. The cluster is not idiomorph-like. Two of its genes have a weak sexM hit (E about 1e-7); I do not treat that as MAT evidence.

Final coverage: 50 target-phylum genomes remain unscanned (fragmentary assemblies, no genome file, or small N50) and 3 more (*Spraguea* sets) were run but are untestable: 53 genomes with no usable result, and the 1.9 Mb *Synchytrium* assembly (11 genomes in total went through Augustus). The unscanned set is mostly Microsporidia (about 25 after the five run here), 13 Kickxellomycota fragments, and a few Entomophthoromycota, Chytridiomycota and Sanchytriomycota genomes. The lineage conclusions of this note hold for scanned genomes only. A test of the rest needs BFD-style annotation (funannotate) or RNA-based gene models, and assemblies with an N50 below 20 kb are unlikely to help even then.

Validation of the strategy: step 2 (structural verification) was applied in 5 lineages. In Taphrinomycotina it removed the Protomyces leads, in Allomyces it found no excess against a random baseline, in *Globomyces* it explained an apparent difference as a tandem duplication, and in the other new lineages (Rozellomycota, Monoblepharomycota, Microsporidia, Calcarisporiellomycota) it could not test the leads across genomes, because each lineage has one genome or the arrangement looks like a tandem family. For those lineages the claims stay "no evidence", not "absent".

### Why the last 50 genomes are not scanned (checked 2026-11-14)
`unscanned_assembly_stats.tsv`: of the 50 genomes not run, the median N50 is 5.4 kb (maximum 41 kb) and the median assembly has 938 contigs; one has no genome file. Only three have an N50 above 15 kb: a 0.5 Mb *Vairimorpha* sp. fragment (12 contigs), the *Pleurotus ostreatus* x *Saccharomyces cerevisiae* fusant (not a lineage without a defined MAT) and *Metchnikovella incurvata* (4.6 Mb, N50 17.6 kb). A flank-first scan needs about 11 genes on one contig (a window of 5 on each side), which is a contig of 20 kb or more in these gene-poor genomes. Annotating these assemblies (funannotate or otherwise) would give gene models but not the contiguity the method needs. This is a limit of the assemblies, so the coverage statement stays: no usable result for these 50 genomes.

## Genome-level matPi search: the PM-like pair is present in 6 of 17 genomes (added 2026-11-14; corrects statements above)
The BFD gene models do not include the homeodomain gene of this locus (the record of *T. deformans* says the same). The curated *T. deformans* matPi is a 181 aa partial protein from the old assembly, and its homeodomain region has no gene model in the BFD annotation either. So an annotation-based search misses it. I searched the assemblies directly.

Scripts: `matpi_tblastn.py` (tblastn of the 4 curated matPi proteins against each assembly), `region_pfam_check.py` (six-frame translation of a region and Pfam), `matpi_hmg_adjacency.py` (distance of every HMG-box gene to the homeodomain-like regions, with a permutation null). Jobs: tblastn and region checks run on the HPCC 2026-11-14; tests in `tests/test_region_pfam_check.py` (4 tests). Tables: `matpi_genome_search_hits.tsv`, `matpi_region_pfam.tsv`, `matpi_hmg_adjacency_pairs.tsv`, `matpi_hmg_adjacency_summary.tsv`, `matpi_hmg_adjacency_null.tsv`. 17 genomes searched (4 *T. deformans*, other *Taphrina* 5, *Protomyces* 7, Savitreella, Saitoella; *T. deformans* JCM22205 has no genome file in the library).

- Positive control: in all three *T. deformans* genomes with a genome file the curated matPi matches an unannotated region at 100% identity (E 1e-119). In A2 it lies on the same scaffold as the matMc gene, 40,370 bp away (RHGF01000006.1: matPi 51,688 to 52,236, matMc 92,606 to 93,298). In CBS 355.35 and PYCC 5710 the two sit on different contigs.
- The top-ranked matPi-like region lies on the same contig as an HMG-box gene in:

| genome | HMG gene | region (identity, E) | distance | Pfam in the region (six-frame ORF) |
|---|---|---|---|---|
| *Protomyces inouyei* | F7C0A940_005637 | 45.2%, 2e-18 | 7,165 bp | Homeodomain, 167 aa, E 1.1e-10 |
| *P. lactucae-debilis* NRRL YB-4353 | F21A9957_003986 | 45.2%, 2e-18 | 7,164 bp | Homeodomain, 151 aa, E 1e-10 |
| *P. lactucae-debilis* 12-1054 | F003EF48_006292 | 45.2%, 2e-18 | 7,152 bp | Homeodomain, 151 aa, E 1e-10 |
| *Taphrina confusa* | F57803F5_005977 | 46.5%, 2e-46 | 14,057 bp | Homeodomain, 168 aa, E 1.1e-8 |
| *T. communis* | F8DCC00A_005608 | 52.4%, 8e-54 | 18,060 bp | Homeodomain, 184 aa, E 1.5e-9 (checked in the first pass, not as the adjacent pair) |
| *T. betulina* | FC6B0312_004634 | 43.0%, 2e-26 | 10,297 bp | no Pfam ORF of 40 aa or more in the region: weaker |
| *T. deformans* A2 (control) | F891DC7E_001655 | 100%, 1e-119 | 40,370 bp | Homeodomain, 189 aa, E 1.5e-9 |
Negative control, a 549 bp stretch 20 kb from the A2 hit: no Pfam ORF. *Savitreella* has a pair at 39 kb (not within 20 kb).
- Permutation null (`matpi_hmg_adjacency_null.tsv`): an HMG gene within 20 kb of a matPi-like region in 6 of 17 genomes; random placement of the same number of HMG genes at random gene positions gives a mean of 0.71 genomes (maximum 5 in 2,000 draws), p = 0.0005 (the floor with 2,000 draws). The genomes are not independent: there are two *P. lactucae-debilis* assemblies, and the three Protomyces offsets are the same (7.15 to 7.17 kb), which is what a conserved locus would give but also what one lineage would give.
- The 20 kb result includes genomes with no detector call (*T. communis*, *T. betulina*): the detector missed them because only one of the two genes had a model.
- Not adjacent: *P. macrosporus*, *P. arabidopsidicola*, *P. gravidus*, *P. pachydermus*, *Taphrina* sp. SM11, *T. wiesneri*, *Saitoella*, and two *T. deformans* assemblies (contigs may be broken at the locus).

What this means, with limits. A homeodomain gene is present in these genomes whether or not an HMG gene is next to it (for example *P. gravidus*, *T. communis* have unannotated homeodomain ORFs), so the homeodomain gene alone is not MAT evidence. The evidence is the conserved physical pairing of an HMG-box gene and a homeodomain gene, which is the arrangement of the curated *T. deformans* locus (matMc plus matPi) and of the Pneumocystis PM locus. The HMG genes in Protomyces are only 44 to 49% identical to the curated matMc and could be any HMG gene; the pairing, not the HMG identity, carries the evidence. I did not test that the HMG genes are orthologous across genomes, I did not look at expression, and the mating-type role of this locus in these species is not shown. The candidate set, for a curator: *P. inouyei*, *P. lactucae-debilis* (two assemblies), *T. confusa*, *T. communis*, and *T. betulina* (weaker).
Next test: proper gene models of both genes in these genomes (for example exonerate of the curated matPi/matMc against the 100 kb region, or the Pfam ORFs of the table), then a check that the homeodomain genes of the locus are orthologous across the genomes and differ between the two *P. lactucae-debilis* strains if they are alleles.
