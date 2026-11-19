# Weak-fragment cap: rank the genes of the called idiomorph (#341)

Date 2026-11-19. Curator decision 2026-11-12 (issue #341). The cap stays opt-in (`--weak-fragment-cap-policy on`), in no preset.

## Change
`weak_fragment_cap.evaluate` ranks only the live core genes whose name starts with the called idiomorph plus `-` (`MAT1-1-` or `MAT1-2-`), by polish status and then bitscore, as before. A call with no live gene of its own idiomorph is not capped (nothing is lowered without evidence about the called side). The note gains `ranked_idiomorph`. Code: `src/MATPredict/detect/weak_fragment_cap.py`; tests: `tests/detect/test_weak_fragment_cap.py`.

## Re-measurement (job 32219550, commit dc3ae31; script `scripts/weak_fragment_cap_remeasure.py`, test `tests/test_weak_fragment_cap_remeasure.py`)
Both rules applied to the frozen v0.9.0 reports (`--preset v0.9.1`, cap off in the reports, so the tier is the final one). Tables: `results/2026-11-19_weak_fragment_cap_341/{ascomycota_v090,fungidb_v090}/{summary,changes}.tsv` on the HPCC.

| Run | Rule | Loci with a called idiomorph | Medium/high | Would change tier | from medium | from high |
|---|---|---|---|---|---|---|
| Ascomycota (reports_all.tar.zst, 2026-11-09_v090_ascomycota) | old | 11,980 | 11,906 | 66 | 20 | 46 |
| same | new | 11,980 | 11,906 | 66 | 20 | 46 |
| FungiDB (313 genomes, 2026-11-09_fdb_v090) | old | 193 | 193 | 3 | 1 | 2 |
| same | new | 193 | 193 | 0 | 0 | 0 |

- Ascomycota: the same 66 loci are matched under both rules (the v0.8.1 run gave 68; v0.9.0 differs by the preset's caps and the SEG-off reruns). One locus changes the deciding gene, none changes tier: GCA_022385605.1 XylFL0594 (MAT1-2, medium): old MAT1-1-3 polished_single 30.4 bits, new MAT1-2-1 unpolished 40.4 bits; capped either way.
- FungiDB: the three loci the old rule capped (*A. terreus* NIH2624 high, *C. geophilum* 1.58 medium, *S. schenckii* 1099-18 high; each called MAT1-1 with a thin MAT1-2-1 fragment, 37.7 to 45.4 bits, as the deciding gene) are kept at their tier. Each has an annotated MAT1-1 gene (note 2026-11-12_290-annotation-by-tier). Every other FungiDB locus was unaffected.
- No locus is newly capped by the redefinition on either run.

## Limits
The FungiDB annotation was not re-run against the new rule: since no tier changes in that run, there is nothing new to compare. The 66 Ascomycota loci that the cap would lower have not been reviewed against annotation; that is a separate step before the cap enters any preset.
