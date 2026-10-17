# Optional-gene policy (`--idiomorph-policy core_first`): campaign comparison

PR #89. Date 2026-10-17. Code: `idiomorph.py` (`resolve_idiomorph_overlaps`, `idiomorph_candidates`, `assign_idiomorph`), threaded through `run_pipeline`, the relaxed pass and the CLI; recorded in `run.parameters.idiomorph_policy`.

## What the policy does
`legacy` (default): unchanged. `core_first`: (a) a mutually exclusive pair in which either gene is `optional` is not resolved, so an accessory gene cannot supersede a core hit; (b) the vote ignores optional genes whenever any core informative gene has a live hit anywhere in the locus (an optional-only locus still votes). Reading (b) per idiomorph, as first worded, would still let MAT1-1-3 outvote MAT1-2-1, so the test is locus-wide; this is also what the 2026-10-16 replay measured.

Other callers, decided explicitly: `resolve_idiomorph_by_models`, `resolve_idiomorph_by_hmm` and the HMM classifier only run for families with `model_idiomorph_alternatives` / `idiomorph_classifier` (Mucoromycota; its genes are not optional), so they are untouched. `evidenced_idiomorphs` (rescue scope, a conservative presence union) is untouched. Tiering reads the call, not the vote. Basidiomycota has 10 optional genes; the policy would apply there if asked, and was not measured.

## Design of the comparison
Code = origin/main at the time of branching (before PR #95, so no partial-core rule), only `--idiomorph-policy` differing; legacy and core_first ran as two parallel processes per genome (same genome, same taxid, MATPREDICT_HTML=0). 8,938 genomes the fix could change by construction plus 500 random controls outside that set (seed 1, drawn from non-set genomes with an Ascomycota:MAT locus in v0.6.0). 28 genomes had no taxid in the assembly summaries: 25 fix, 3 control dropped. Run: 8,913 fix + 497 control; 18,820 runs, all rc=0. Output in /bigdata/stajichlab/jstajich/projects/optgene_runs.

## Results
Fix set (8,913 genomes): loci 9,321 (legacy) vs 9,468 (core_first). 473 genomes (5.3%) change label or gain/lose a locus. Locus-level label changes:

| change | loci |
|---|---|
| absent -> MAT1-2 (locus gained) | 202 |
| MAT1-1 -> MAT1-2 | 177 |
| MAT1-2 -> absent (locus lost) | 71 |
| MAT1-1 -> undetermined | 46 |
| absent -> MAT1-1 | 18 |
| MAT1-1 -> absent | 2 |

Confidence high->medium 1,426, medium->high 165, low->medium 8. Locus class mat_locus -> homothallic_candidate 842, partial_locus -> mat_locus 96 (both cores now live at the same locus).

By order (genomes with a label/loci change): Hypocreales 110, Diaporthales 105, Mycosphaerellales 70, Sordariales 39, Eurotiales 34, Onygenales 15, Lipomycetales 11, Dipodascales 10. By genus: Fusarium 47 (33 MAT1-1->MAT1-2 loci), Cryphonectria 45 (all MAT1-1 -> undetermined), Diaporthe 38 (37 MAT1-1->MAT1-2), Neurospora 24, Aspergillus 22 (13 lost, 9 gained), Cordyceps 17, Cercospora 13 and Fulvia 11 (gains), Microsporum 11 and Lipomyces 9 (MAT1-1->MAT1-2).

Control (497): NOT identical. 11 genomes (2.2%) change, all but one a locus gained (Pyricularia oryzae and relatives, MAT1-2-1 + MAT1-1-3 cross-hit: legacy superseded MAT1-2-1, the one remaining gene failed the evidence floor, the locus was dropped). So the by-construction selection missed loci that did not exist in the v0.6.0 reports; genomes with no locus in v0.6.0 were not run at all and may gain loci the same way. The size of that is unmeasured.

### Agreement with independent evidence
| evidence | legacy | core_first |
|---|---|---|
| BioSample truth, 203 genomes (144 in run set; others v0.6.0 label, unchanged by construction): agree | 178 | 150 |
|   wrong / undetermined / mixed / no locus | 21 / 0 / 0 / 4 | 23 / 25 / 1 / 4 |
| Alpha/HMG profile calls (PR #86 scored.tsv; loci present in both runs) | 444 of 607 | 522 of 604 |
| Parastagonospora, 109 genomes in run set, with records | MAT1-2 91 (88 with both core genes), MAT1-1 18 | identical |

Truth set: no genome improves; 28 get worse, in three species: Cryphonectria parasitica 25 genomes (stated MAT1-1; legacy MAT1-1, core_first undetermined), Juglanconis oblonga 2 (MAT1-1 -> MAT1-2, wrong), Botrytis cinerea Sl8 (MAT1-1 -> mixed; the known 16%-coverage MAT1-2-1 relic, PR #95). The 25 C. parasitica genomes are one clonal species, so it is closer to 3 independent regressions than 28.

Mechanism of the regressions: the strong evidence at those loci is on an optional gene (MAT1-1-2 142 bits in C. parasitica; MAT1-1-1/-2/-3 up to 149 in Juglanconis) while the core MAT1-1-1 hit is weak (55.8 and 39.7 bits); ignoring optional genes leaves a weak core MAT1-1 tied with or below a weak MAT1-2-1 (HMG-like) hit. The rule discards real signal when the core reference is diverged.

Parastagonospora: the v0.6.0 failure (87 loci, MAT1-1) is already gone in legacy on current main (MAT1-2 at 91 of 109 with records), so core_first adds nothing there; I did not trace which post-v0.6.0 change fixed it. Of 9,812 v0.6.0 Ascomycota:MAT loci in the run set, 7,875 have a locus at the same contig/start in legacy-on-main and 7,789 of those carry the same label; unmatched loci differ in coordinates.

Runtime: 452.3 vs 456.2 CPU-hours (+0.9%) summed over 9,410 genomes, median 165 vs 167 s.

## Reading
Profile agreement improves (73% -> 86%) but is not independent evidence: the profiles score MAT1-1-1/MAT1-2-1 entries only, so they share core_first's blind spot to optional genes. The only independent check, BioSample, moves from 178 to 150 agreeing with no gains. 473 of 8,913 genomes change, concentrated in Fusarium, Cryphonectria, Diaporthe, Neurospora, Aspergillus, Fulvia, Cercospora, Cordyceps; for most of the MAT1-1 -> MAT1-2 changes (Fusarium, Diaporthe, Microsporum, Lipomyces, Neurospora) there is no independent truth in the set, so whether they are corrections is unverified. The locus gains (220, Pyricularia, Fulvia, Cercospora, Cordyceps, Hortaea) look like real MAT1-2 loci rescued from supersession but are not truth-checked.

Candidate refinement (not implemented): apply only part (a), or let optional genes vote when the core hit is weak/low-coverage (e.g. PR #95's coverage rule) rather than never.

## Assumptions and not verified
- Campaign code is main before PR #95; the branch was since merged with main and its default behaviour checked identical to main on 2 genomes (reports equal excluding `run`).
- Truth agreement is per genome: all Ascomycota:MAT loci share one label; otherwise "mixed". Genomes with several loci can differ from per-locus scoring.
- Profile agreement uses v0.6.0 scored.tsv profile calls keyed by (genome, contig, start); they were not rescored on the new reports.
- Parastagonospora "expected" is not defined by a record; I report call and presence of both core genes.
- HTML and diagnostics off; the 28 taxid-less genomes were not run; legacy and core_first ran concurrently so wall times are comparable but not isolated.

# Follow-up: `resolution_only` (part (a) of core_first only)

Curator's decision after the first comparison: `core_first` is not merged as is. Third policy value `--idiomorph-policy resolution_only`: a pair in which either gene is `optional` is not resolved; the vote is exactly the legacy vote (optional genes keep voting). Unit tests in tests/detect/test_idiomorph_policy.py (Parastagonospora pattern: both stay live and the vote decides; Cryphonectria pattern: strong optional MAT1-1-2 with weak core hits still calls MAT1-1 under legacy and resolution_only, but undetermined under core_first). tests/detect: 966 passed locally, including the four-combination partial-core test.

## Run
Code = branch merged with current main (includes PR #95; `--partial-core-policy` off), legacy vs resolution_only as parallel processes per genome. Sets: (1) the 473 genomes that changed under core_first (470 run), (2) the 497 controls, (3) the 203 BioSample truth genomes, (4) 300 random genomes (seed 1) with `uncalled` status and lineage/phylum_fallback routing in v0.6.0, i.e. no locus (298 run). 1,429 of 1,430 genomes complete in both policies (one task was still running at harvest; 2 set-4 genomes and the dropped taxid-less ones excluded). All rc=0.

| set | genomes | loci legacy / res_only | genomes with label or loci change |
|---|---|---|---|
| 1 changed under core_first | 470 | 601 / 746 | 272 |
| 2 controls | 497 | 504 / 514 | 10 (all a locus gained) |
| 3 BioSample truth | 203 | 200 / 200 | 0 (20 with confidence/class changes) |
| 4 no locus in v0.6.0 | 298 | 14 / 78 | 65 (64 gain a MAT1-2 locus) |

Truth (203, stated MAT1-1/MAT1-2): legacy 178 agree / 21 wrong / 4 no locus; resolution_only identical, genome by genome. Cryphonectria parasitica (25), Juglanconis oblonga (2) and Botrytis Sl8 are unchanged from legacy. Not worse anywhere; also not better anywhere.

Label/loci changes over all sets (labels at matched coordinates): loci gained 292 (2 more genomes already had loci), lost 73, MAT1-1 -> MAT1-2 only 20 in set 1 (core_first: 177 in the full fix set), MAT1-1 -> undetermined 0. Gained loci by genus: Fulvia 22, Fusarium 20, Cordyceps 15, Cercospora 13, Komagataella 13, Pyricularia 10, Hortaea 10, Friedmanniomyces 9, Sporothrix 8, Aspergillus 8. Of 292 gained, 267 carry MAT1-2-1 + MAT1-1-3 (the Pyricularia pattern); 207 of the 227 gains in genomes also run under core_first were gained there too, so resolution_only keeps the useful part of core_first. Controls: Pyricularia 7, Ceratocystis, Lasiosphaeris, Phyllosticta (same 10 as before).

Caveat on the gains: 198 of the 292 gained loci have class `idiomorph_gene_only` (MAT1-2-1 + MAT1-1-3 hits and no flanking gene), medium confidence, 44 are high/mat_locus. Pyricularia oryzae MAT1-2 is plausible; 13 Komagataella (Pichiales) and Xylaria/Geotrichum/Saccharomycopsis gains in set 4 are not truth-checked and could be HMG-box proteins cross-hitting MAT1-2-1 rather than MAT loci. Set 4: 65 of 298 genomes (22%) gain one locus: the policy matters most among genomes with no locus, and that is also where false positives would live.

Lost loci (73): 68 are partial_locus MAT1-2-1 + flank loci (Aspergillus 13, Neurospora, Neosartorya...); 72 of 73 are in genomes that still have a locus in resolution_only (e.g. on another contig), only 2 coincide with a gained locus within 30 kb. I did not trace the mechanism (candidates: polish-cap rank or the better-supported cluster displacing a partial one).

Confidence: set 1 high->medium 12, medium->high 14, low->medium 7; mat_locus -> homothallic_candidate 55, partial_locus -> mat_locus 50. Truth set: high->medium 17, medium->high 3, mat_locus -> homothallic_candidate 16 (no label change).

Parastagonospora: no change (18 in these sets, 17 MAT1-1 with both core genes; the 109-genome run earlier showed legacy already right on main). Alpha/HMG profile agreement (v0.6.0 scored.tsv, not independent: profiles score MAT1-1-1/MAT1-2-1 only): legacy 57/154, resolution_only 55/151 on loci present in both; these sets are enriched for contested loci, so the rates are not comparable with the earlier 73% / 86%.

Unverified: whether the gained loci are real; the cause of the 73 lost; one genome not finished; resolution_only on Basidiomycota; interaction with `--partial-core-policy coverage_profile` beyond the unit test.
