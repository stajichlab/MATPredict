# Ascomycota, v0.7.0 definitive run (`cleanrun-2026-10c`): judgement

The definitive run's Ascomycota calls are the first clean run's calls. Of 19,343 reports, 19,338 are identical to `cleanrun-2026-10b`
in calls and gene models; 5 differ in a polished gene model only (calls unchanged): GCA_030564055.1, GCA_030569375.1,
GCA_057550925.1, GCA_060041085.1 and GCF_027921745.1 (the same gene is placed on a different one of two equal-scoring copies;
same code run twice gives different models, issue #161). The counts against v0.6.0 are those of the first run's judgement below: 150 calls
lost, 472 new, 588 confidence down, 387 up, 219 idiomorph flips, 759 `label_unsupported` loci. The tables are in
`results/2026-10-26_v070_ascomycota/`. The judgement below was written for `cleanrun-2026-10b` on 2026-10-23 and is kept as written; its
verdict applies to the definitive run through the identity above. One addition: the defective-gene flag (#152, report-only) was not in 10b; in
this run it marks 5,160 of 20,985 loci, mostly low-identity flank-gene models in distant orders (pilot: 7.5% of polished entries), and its
table presentation is a curator decision.

---

# Ascomycota clean run (cleanrun-2026-10b) vs v0.6.0: phase-5 judgement (2026-10-23)

**Inputs.**
- Clean run: `results/2026-10-23_ascomycota_cleanrun`, 9 waves, 19,415 listed genomes, 19,343 reports. Tag `cleanrun-2026-10b`, commit 56e11cd, the five `DETECT_ARGS` flags as in `MANIFEST.txt`.
- Baseline: the v0.6.0 campaign (`.claude/worktrees/run-56e11cd/results/2026-10-03_ascomycota_v060`). All 19,380 of its per-genome reports were extracted from `reports_all.tar.zst`.

**Where the diagnostics are.** Everything is on the HPCC under `results/2026-10-23_ascomycota_diag/`:
- `all_genomes_final.tsv`: one row per listed genome, with its category and old and new calls.
- `changes/<category>.tsv`: the per-genome change lists.
- `by_class_categories.tsv`, `by_order_categories.tsv`.
- `attr.tsv`: report-evidence signatures for every changed genome.
- `confidence_DOWN_mech.tsv`, `conf_up_mech.tsv`, `mat_accounting.txt`, `truth_check.txt`, `high_split.tsv`, `label_unsupported.tsv`.
- The re-runs: `rerun_sample*.tsv`, `runs_<config>/`, `rerun_attribution*.tsv`, `rerun_summary.txt`.
- `scripts/`, plus `old/` (the v0.6.0 reports) and `old.jsonl`/`new.jsonl` (compact summaries).

**What was not touched.** No repo, frozen worktree, wave output or Basidiomycota job was edited or touched.

**Compute.** Re-runs were 14 `highclock` array tasks of 16 CPUs each (jobs 30344986 and 30354356), about 1,590 detect runs in all. Parsing was done on the login node with `nice`. The first attempt to load everything at once was killed by the login-node memory limit and was rewritten to stream.

## Verdict

**The campaign passes phase 5 with a short curator list. No change goes against independent evidence.**

Evidence checks:
- **BioSample truth (296 genomes).** 16 move from wrong to right: 15 *Botrytis cinerea* and 1 *B. fabae*, all MAT1-2, fixed by the partial-core rule. One *Sclerotinia minor* stated "homothallic" moves from a single label to `undetermined`. No truth genome moves from right to wrong.
- **Dothideomycetes box-ORF scan** (the 2026-10-18 audit, `gw_class`):
  - 87/87 *Parastagonospora* MAT1-1 -> MAT1-2 flips now agree with "no alpha-box ORF in the genome" (before, all 87 disagreed).
  - 39 of 40 scanned NEW calls agree.
  - One changed genome goes from agree to disagree (below).

Attribution:
- **Re-runs.** 546 sampled genome-category rows: 492 in the main batch (482 genomes; every large category has at least 40) plus 54 supplementary genomes in the largest clade groups. Of these, 16 were reproduction controls: v0.6.0 code + v0.6.0 db reproduced the v0.6.0 campaign in 16/16, and the frozen code + flags reproduced the clean run in 16/16, so the diffs are not noise.
- **Every category has a cause.**
  - Flag-driven: the homothallic policy, the partial-core rule and the scoped roster, each confirmed by toggling its flag.
  - Database-driven: most of the rest. In these, v0.6.0 code + the new db reproduces the change and the new code + the old db does not.
  - Code-driven: a small remainder, namely the genetic-code-26 fix and the subphylum guard (the guard needs both the new code and the new db).

**The MAT1-1 drop is fully explained** (section 2).

**Three systematic tier artefacts account for most of the 588 confidence-DOWN calls.** They are reporting issues, not wrong calls:
- the #137 margin fallback: 380, of which *Macrophomina* 228 and *Pyricularia* 98;
- a distant unpolished fragment that is no longer ignored by the tier: 93 *Pyricularia*;
- the *Aureobasidium* vote collapse: 30.

## 1. Category table (whole campaign, all 9 waves)

Rules:
- Categories are assigned per genome, using the pilot's precedence and its primary locus (the first sorted locus).
- A coordinate-only shift under 1 kb counts as identical, as in the pilot. That covers 2,925 genomes, median shift 6 bp.
- "secondary locus changed" is new here: the primary locus is identical, but a second locus changed idiomorph or confidence.
- Pilot rates are for the 461 Ascomycota pilot genomes. The pilot list was stratified and enriched for problem clades, so its rates are higher by construction.

| category | n | % | pilot Asco % | cause (evidence) | verdict |
|---|---|---|---|---|---|
| identical call | 13,171 | 67.8 | | | |
| uncalled both | 1,601 | 8.2 | | | |
| identical + uncalled | 14,772 | 76.1 | 58.6 | | |
| no new report | 72 | 0.4 | 0 | All 72 are on the BFD suppress list. 44 had v0.6.0 calls (93 loci): *S. cerevisiae* 38, *Nakaseomyces glabratus* 2, *Blumeria graminis* 2, *F. graminearum* 1, 1 unassigned. For example, GCA_003277255.1 is "flagged low quality at NCBI (user decision 2026-10-05)". 16 were uncalled and 12 had no v0.6.0 report. | expected (curation decision, not a detection change) |
| NEW call | 472 | 2.4 | 14.5 | <ul><li>160 *Komagataella* `MATkp` (new family). All MATa only; the MATalpha cassette is withheld by `modelled_gene_bar` in all 160 (see 6).</li><li>111 *Pleosporaceae* scoped roster: *Bipolaris* 56, *Curvularia* 23, *Alternaria* 14, *Exserohilum* 10, *Stemphylium* 7.</li><li>51 *Geotrichum* `MATdc`, 28 *Yarrowia* `MATyl`.</li><li>16 genetic-code-26 genomes (Alaninales etc.) with no v0.6.0 report.</li><li>About 105 new-record MAT calls: *Seiridium* 16, *Ceratocystis* 15, *Zymoseptoria* 12, *Corynespora* 8, *Pestalotiopsis* 8, *Apiospora* 6, others.</li></ul>Re-runs (48): db 23, scoped roster 6, code 4 (*Zymoseptoria* Zt117, *Aspergillus*, *Brettanomyces*, *Exophiala*), code-26 3 (old code exits on `--geneticcode 26`), new code + new db 3, partial 7, either alone 2. | expected / improvement. Box-ORF agrees in 39/40 scanned Dothideomycetes. |
| LOST call | 150 | 0.8 | 2.0 | 113 subphylum guard: the old family is not attempted. Lipomycetes 51, Dipodascomycetes 36, Pichiomycetes 12, Saccharomycetes 7, Trigonopsidomycetes 5, Orbiliomycetes 2. 37 attempted: *Eutypa* 12, *Phialophora* 5, *Neopestalotiopsis* 3, *Schizothecium* 2, *Metschnikowia* 2, and 5 Hypocreales (*Cordyceps* 2, *Ophiocordyceps*, *Simplicillium*, *Clonostachys*). Re-runs (40): new code + new db 15 (guard), db 17, partial 8 (guard + db). | guard: **expected** (plan §6.4 expected >= 52 Lipomycetales + 20 Dipodascales; observed 53 + 36). Weak attempted losses: expected. **5 Hypocreales: needs attention.** |
| idiomorph -> undetermined | 626 | 3.2 | 6.5 | <ul><li>**606 homothallic policy.** In 597 `vote_idiomorph` equals the old call (*Fusarium* 430, all *F. graminearum* species complex: *F. asiaticum* 253, *F. graminearum* 128, ...; *Pyrenophora tritici-repentis* 62, *Sclerotinia* 24, *Ceratocystis* 16, *Neurospora* 15, *Didymella* 13, *Aspergillus* 10, *Stemphylium* 5). In 9 the vote differs from the old call because new records moved it (*Xanthoria* 2, *Rusavskia*, *Botrytis*, *Botryotinia*, *Stromatinia*, *Calogaya*, *Leptosphaerulina*, *Polytolypa*).</li><li>19 partial-core `undetermined`, no homothallic note: MTL yeasts *Schwanniomyces*/*Debaryomyces* 13 (already `homothallic_candidate` in v0.6.0), *Sclerotinia* 3, *Calonectria* 2.</li><li>1 *Aureobasidium*: tie, db.</li></ul>Re-runs (69): `nohomo` reverts 36/40 main-group genomes. The other 4 revert the idiomorph but differ in class or confidence from concurrent db changes. Partial-core 16. | **expected** (biology: FGSC, Ptr, *S. sclerotiorum*, *Ceratocystis* are homothallic or unisexual) |
| idiomorph flip | 219 | 1.1 | 4.1 | <ul><li>137 MAT1-1 -> MAT1-2: *Parastagonospora* 87 (db, new `13684` records), partial-core 40 (*Botrytis* 30, *Ophiostoma* 5, *Sclerotinia*, *Botryotinia*, *Leptographium*, *Intubia*, *Sporothrix*), *Xylaria* 3, *Ustulina* 2, ...</li><li>40 MAT1-2 -> MAT1-1: *Aureobasidium* 17, *Friedmanniomyces* 6, *Xylaria* 3, *Thielaviopsis* 2, ...</li><li>33 yeast relabels into the new scope-only families (`MATdc` A/B, `MATog`/`MATkp` MATa, `MATyl` B).</li><li>6 undetermined -> MAT1-1 (*Phyllosticta capitalensis* 5).</li></ul>Re-runs (60): db 46, partial-core 9, partial 5. | *Botrytis*/partial-core: **improvement** (BioSample 16/16). *Parastagonospora*: **improvement** (box-ORF 87/87). Yeast relabels: expected. ***Aureobasidium*, *Friedmanniomyces*, *Phyllosticta capitalensis*, *Xylaria*: needs attention** (label unsupported, below). |
| confidence up | 387 | 2.0 | 2.4 | 114 GAP1/BGL1 newly polished (*Alternaria* 106, *Drechslera* 3; scoped roster). 84 an unpolished opposite gene now ignored or gone (*Ascochyta*, *Corynespora*, *Plenodomus*, ...). 71 a core gene newly polished (*Diaporthe*, *Trichoderma*, *Pseudogymnoascus*, ...). 62 margin now >= 5. 41 `partial_locus` -> `mat_locus`. Re-runs (40): db 36, scoped roster 2, code 1, either 1. | expected |
| confidence DOWN | 588 | 3.0 | 1.1 | <ul><li>**380 #137 single-candidate identity-margin fallback < 5:** *Macrophomina* 228 (margin 5.18 -> 4.96 from a 0.2-point shift in MAT1-2-1 identity), *Pyricularia* 98, *Trichoderma* 15, *Fusarium* 14, *Alternaria* 6, ...</li><li>**93 *Pyricularia*:** a distant (36 kb) unpolished MAT1-1-2 fragment, plus the MAT1-1-3 overlap loser, is no longer in `confidence_ignored_genes`. The margin *rises* (33 -> 36.6 bits), yet the tier drops to medium and the span grows 14 -> 51 kb.</li><li>55 a polished gene becomes unpolished: MTLA1 in CTG yeasts (*Metschnikowia* 18, *Scheffersomyces* 7, *Candida* 5, ...).</li><li>32 bit-margin collapse < 5: *Aureobasidium* 30.</li><li>16 class also changed; 5 other.</li></ul>Re-runs (48 + 23 supplementary): db-driven 47/48. The one exception is a homothallic conversion of the other locus (*F. graminearum*), a primary-locus artefact. | artefacts of the tier (#137, #129-type distant fragment): **needs attention as a reporting rule, not as wrong calls**. *P. oryzae* and *M. phaseolina* are heterothallic and the MAT1-2 label is unchanged. |
| class changed | 214 | 1.1 | 0.2 | <ul><li>196 *Aspergillus* + 5 *Neosartorya*: `homothallic_candidate` -> `mat_locus`, MAT1-1 high. The `af293` record correction makes the in-locus MAT1-2-1 a 363-bp, 37.6%-coverage fragment, as in the pilot. 180 are *A. fumigatus*.</li><li>13 `partial_locus` -> `mat_locus` (*Podospora* 5, *Schizothecium* 3, ...).</li></ul>Re-runs (40): db 40/40. | *A. fumigatus*: **improvement** (heterothallic). **11 genomes need attention:** *Neosartorya* sp. x5, *A. hiratsukae* x2, *A. spinosus*, *A. latus*, *A. felis*, plus 3 *A. fumigatus*. These carry a separate full MAT1-2 locus, so the genome holds both idiomorphs (homothallic species or a hybrid such as *A. latus*). |
| span changed >1 kb | 1,300 | 6.7 | 8.0 | 453 same genes and `core_span` equals the old extent (outer-cluster artefact). 331 same genes, `core_span` differs. 155 GAP1 added (scoped roster). 114 MAT1-1-1, 90 MAT1-1-2, 40 MAT1-2-1 and 65 MAT1-2-7 added (new records; distant fragments). By genus: *Colletotrichum* 149, *Alternaria* 138, *Pyricularia* 116, *Fusarium* 66. Re-runs (40): db 22, new code + new db 8, scoped roster 3, other 7. | cosmetic / expected. Maps should draw `core_span`. |
| locus count changed | 306 | 1.6 | 3.5 | Spurious `PM` loci removed in *Pichia* (76; subphylum guard / scope). *Ogataea* -> `MATog` (24). *Komagataella* -> `MATkp` (15). `MATyl` duplicates removed in *Petasospora*, *Williopsis*, *Ogataea* (17). *Aspergillus* 1 -> 2 (13; a second locus). Re-runs (40): db 20, new code + new db 10, partial 9, code 1. | expected |
| genes_found changed only | 298 | 1.5 | 2.0 | MAT1-2-1 added to MAT1-1 calls (121, *Parastagonospora* 86: an unpolished 58-aa cross-hit, as in the pilot). MAT1-1-3 dropped (72, *Pyricularia* 48). MAT1-2-7 added (32, *Ceratocystis*). 19 a gene moved between two loci of the same genome. Re-runs (40): db 32, partial 7, scoped roster 1. | cosmetic |
| secondary locus changed | 11 | 0.1 | 0 | *Friedmanniomyces* x3 (MAT1-2 -> medium), *Aureobasidium*, *Agyrium*, *Hortaea*, *Millerozyma*/*Saccharomyces* (homothallic / partial-core undetermined), *Williopsis*, *Ambrosiozyma*. Re-runs: db 7, new code + new db 4. | expected |

**Whole-campaign totals.** Called genomes rise from 17,314 to 17,592, a gain of 278. That is 472 NEW, minus 150 LOST, minus 44 suppressed genomes that were called in v0.6.0. The 472 NEW include 16 genomes that failed in v0.6.0 (genetic code 26).

**Correction made during this pass.** The first `cmp.py` left the v0.6.0 calls of no-report genomes blank. It is fixed, and `all_genomes_final.tsv` and `changes/no_new_report.tsv` have been rewritten (`scripts/fix_norep.py`). The locus accounting in section 2 was not affected: there the 44 suppressed genomes fall under "locus lost with the genome".

**By class** (`by_class_categories.tsv`; share of genomes changed):
- Pleosporales are the most-changed order, at 78.7% of genomes changed.
- Dothideomycetes 48.2%: *Parastagonospora*, *Macrophomina*, Pleosporaceae.
- Dipodascomycetes 39.2% (guard, `MATdc`/`MATyl`).
- Sordariomycetes 28.8% (*Fusarium* homothallic, *Pyricularia* tier).
- Leotiomycetes 28.0% (*Botrytis*, *Sclerotinia*).
- Pichiomycetes 17.9% (*Komagataella*, *Pichia* `PM`).
- Eurotiomycetes 12.6% (*A. fumigatus* class).
- Saccharomycetes 6.2%.
- Pezizomycetes, Schizosaccharomycetes, Taphrinomycetes and Pneumocystomycetes: 0%.

**By order** (`by_order_categories.tsv`), largest shifts:

| order | n | changed % | main driver |
|---|---|---|---|
| Lipomycetales | 59 | 94.9 | 53 LOST (guard) |
| Pleosporales | 1,012 | 78.7 | NEW 121, up 203, flips 89, undetermined 91, span 186 |
| Pichiales | 579 | 76.0 | NEW 161, locus count 183 |
| Alaninales | 20 | 75.0 | 15 NEW (code 26) |
| Magnaporthales | 544 | 67.6 | DOWN 191, span 116 |
| Microascales | 164 | 65.9 | *Ceratocystis* |
| Phaffomycetales | 188 | 55.9 | |
| Dipodascales | 398 | 39.2 | |
| Botryosphaeriales | 807 | 36.1 | DOWN 229 |
| Hypocreales | 3,535 | 20.7 | undetermined 432 |
| Eurotiales | 2,489 | 13.5 | class 201 |
| Saccharomycetales | 2,307 | 1.0 | |
| Pezizales | 172 | 0 | |

Rates against the pilot: the campaign is less changed in every call category except confidence DOWN and class changed. Those two are dominated by three clades (*Macrophomina*, *Pyricularia*, *A. fumigatus*) that the pilot sampled only lightly.

## 2. The MAT1-1 drop (-619) vs MAT1-2 (-4)

Method: a locus-level match of each old MAT1-1 or MAT1-2 locus to a new locus on the same contig with an overlapping span. The result is `mat_accounting.txt`, which nets -618 and -2. The 1-2 locus differences come from overlap matching.

| | MAT1-1 | MAT1-2 |
|---|---|---|
| -> `undetermined`, homothallic, `vote_idiomorph` = old call | **-552** (*Fusarium* 437, *Pyrenophora* 62, *Didymella* 13, *Aspergillus* 10, *Neurospora* 9) | -52 (*Sclerotinia* 22, *Ceratocystis* 16, *Neurospora* 6) |
| -> `undetermined`, other (vote != old / partial-core only) | -11 | -4 |
| flipped to the other idiomorph by the partial-core rule | -40 (*Botrytis* 30, *Ophiostoma* 5) | +40 |
| flipped, db-driven | -97 (*Parastagonospora* 87) vs +36 back (*Aureobasidium* 19) | -36 vs +97 |
| relabelled into yeast families (MATalpha/MATa/A/alpha) | -23 | -5 |
| locus lost with the genome (guard, weak losses, suppressed genomes: MAT1-1 3 = *Blumeria* 2 + *F. graminearum* 1; MAT1-2 1) | -52 | -131 |
| locus lost, genome still called | -7 | -59 |
| NEW genome call | +87 | +131 |
| new locus in an already called genome | +28 | +17 |
| from `undetermined` | +13 | 0 |
| **net** | **-618** | **-2** |

**The asymmetry is explained.**
1. The homothallic policy removes 563 MAT1-1 but only 56 MAT1-2. On a fused locus the bitscore vote favours MAT1-1, because MAT1-1-1 plus MAT1-1-2/-3 add up, and the *F. graminearum* complex (430 genomes) is MAT1-1 by that vote. In 597/606 cases the note's vote equals the old call.
2. Flips move about 101 loci net from MAT1-1 to MAT1-2 (*Parastagonospora* 87 db, *Botrytis* etc. 40 partial-core, against *Aureobasidium* and others 36 the other way). That cancels almost all of the MAT1-2 homothallic loss.
3. MAT1-2 loses more to the guard and the yeast losses (131 vs 52). It gains more from NEW calls (131 vs 87: scoped roster, *Seiridium*, *Corynespora*, *Pestalotiopsis*).

## 3. Attribution by re-run (frozen code, flags toggled; v0.6.0 code + new db; new code + v0.6.0 db)

Configurations:
- `off`: frozen code and db, no flags.
- `nohomo` / `nopc` / `noscope`: production flags minus one.
- `v060code_newdb`: worktree run-7c7ed99 code with the clean-run db.
- `newcode_v060db`: frozen code with the 7c7ed99 db.
- `v060` and `prod`: the reproduction controls.

Three `v060code_newdb` runs failed, all `exonerate --geneticcode 26`, which is the known v0.6.0 defect. Verdicts (`rerun_attribution.tsv`, `rerun_summary.txt`):

| category (sampled) | flag-driven | db-driven | code-driven | needs both | other/partial |
|---|---|---|---|---|---|
| confidence DOWN (48) | 1 (homothallic, other locus) | 47 | 0 | 0 | 0 |
| -> undetermined (69) | 36 homothallic + 16 partial-core + 4 combined flags | 9 (vote moved by records, then policy) | 0 | 0 | 4 (`nohomo` reverts the idiomorph only) |
| class changed (40) | 0 | 40 | 0 | 0 | 0 |
| idiomorph flip (60) | 9 partial-core | 46 | 0 | 0 | 5 |
| LOST (40) | 0 | 17 | 0 | 15 (guard: `route()` + `fallback_subphylum`) | 8 (guard + db) |
| NEW (48) | 6 scoped roster | 23 | 4 + 3 code-26 | 3 | 9 |
| confidence up (40) | 2 scoped roster | 36 | 1 | 0 | 1 |
| locus count (40) | 0 | 20 | 1 | 10 | 9 |
| span >1 kb (40) | 3 scoped roster | 22 | 1 | 8 | 6 |
| genes_found only (40) | 1 scoped roster | 32 | 0 | 0 | 7 |
| secondary locus (11) | 0 | 7 | 0 | 4 | 0 |

The 5 kb intron cap (#92) was not the cause of any sampled loss. Only 3 old LOST calls had a gene model over 5 kb (2 *Lipomyces*, 1 *Wickerhamiella*), and all 3 are guard losses.

**Supplementary batch** (54 genomes, job 30354356, `rerun_attribution2.tsv`). The main sampler took only 1 genome from each of the largest single-signature groups, so these were added:

| group | n | result |
|---|---|---|
| *Macrophomina* confidence DOWN | 8 | db 8/8 |
| *Pyricularia* confidence DOWN (4 margin fallback, 4 distant fragment) | 8 | db 8/8 |
| *Aureobasidium* confidence DOWN | 4 | db 4/4 |
| *Metschnikowia* confidence DOWN | 3 | db 3/3 |
| *Parastagonospora* flips | 6 | db 6/6 |
| *Botrytis* flips | 4 | partial-core 4/4 |
| *Aureobasidium* flips | 3 | db 3/3 |
| *Fusarium* -> undetermined | 6 | homothallic 6/6 |
| *Eutypa* LOST | 3 | db 3/3 |
| Hypocreales LOST (*Clonostachys*, *Simplicillium*, *Cordyceps* YNP15) | 3 | db 3/3 |
| *Bipolaris* NEW | 3 | scoped roster 3/3 |
| *Komagataella* NEW | 3 | db-driven (the new `MATkp` record family): v0.6.0 code + new db gives the new call, new code + v0.6.0 db gives the old `MTL`/`MATsc` state |

**Across both batches, 546 sampled rows** (530 change rows plus 16 controls). Every sampled change is attributed to a flag, the db, the code, or new code + new db together (the guard). The "partial" rows are guard + db combinations, or yeast families whose old-db state differs from both ends. None is unexplained.

## 4. *Pleosporaceae* `high` calls: supported vs tie-break (#129 rule with the pilot's `lost_to` amendment)

| | supported | tie-break | grey |
|---|---|---|---|
| `high` idiomorph calls | **296**: *Alternaria* 262 (6 NEW), *Pyrenophora* 28, *Drechslera* 4, *Bipolaris* 1, *Curvularia* 1 | **11**: *Curvularia* 6, *Exserohilum* 5; all 11 are NEW scoped-roster recoveries | **0** |
| `high` `homothallic_candidate` (`undetermined`) | 0 | 5 (*Stemphylium vesicarium* 4, *Decorospora* 1); the split does not apply | 0 |

The pilot had 18 / 6 / 1. *Curvularia* BF9436 was the pilot's grey case; with the amendment it is now tie-break. As predicted, every tie-break `high` call is a new scoped-roster call. All NEW *Bipolaris* (56) are `medium`.

Campaign-wide, over all `high` Ascomycota:MAT calls: 9,163 supported, 29 tie-break and 8 grey. Outside *Pleosporaceae* there are 13 tie-break and 8 grey, mostly *Mycosphaerellales*: *Teratosphaeria* 5, *Cercospora* 4, *Cladosporium* 3, *Pseudocercospora* 2. Each rests on an opposite-idiomorph core fragment (see `high_split.tsv`).

## 5. "Label unsupported" list (#138) (`label_unsupported.tsv`)

The literal rule (genus *Aureobasidium*/*Xylaria*, or `idiomorph_margin` < 20) marks 4,493 of 20,985 loci. But for single-candidate loci `idiomorph_margin` is the identity-point fallback (`pipeline.py` ~l.3292), not bits. Splitting the two:

- **Bit-margin < 20 (>= 2 candidates) or the genus rule: 757 loci.** By genus: *Aureobasidium* 199 (all of its calls), *Zygosaccharomyces* 70, *Botryosphaeria* 50, *Xylaria* 32 (all of its calls), *Friedmanniomyces* 25, *Maudiozyma* 25, *Saccharomycopsis* 23, *Monosporozyma* 21, *Neurospora* 15, *Barnettozyma* 14, *Saccharomyces* 14, *Arxiozyma* 14, *Seiridium* 14, *Pichia* 11, *Kazachstania* 11, *Nigrospora* 9, *Torulaspora* 8, *Phyllosticta* 8, ... By family: MAT 314 + genus rule 231, MATsc 223, MTL 54. By tier: 2,496 + 145 high among all 4,493; for the 757, see the TSV.
- **Single-candidate identity margin < 20: 3,736 loci** (2,308 high), the #137 quantity. By genus: *Beauveria* 364, *Colletotrichum* 350, *Penicillium* 279, *Macrophomina* 239, *Aspergillus* 159, *Pyricularia* 155, *Saccharomyces* 131, *Pichia* 126, *Fusarium* 107, *Candida* 106, *Trichoderma* 103. These are mostly unambiguous single-idiomorph loci: there is no opposite candidate at all, and the "margin" is the identity gap of an HMG cross-hit resolution. **Curator decision needed: whether #138 applies to them.** Recommendation: no. Apply it to bit margins only, plus the genus rule.

## 6. Needs curator attention

1. ***Aureobasidium* (199 loci; 17 flips, 30 confidence DOWN)** and ***Friedmanniomyces*** **(25 loci; 6 flips, 3 secondary).** These are both-gene fragment loci with vote margins of 2-19 bits, and the flips are db-driven. The box-ORF scan says "MAT1-1: HMG ORF at the locus", which fits homothallism. *Aureobasidium* is covered by #138. *Friedmanniomyces endolithicus* is not on the genus list: add it, or rely on the margin rule (it catches 25).
2. ***Phyllosticta capitalensis*** **(5 undetermined/MAT1-2 -> MAT1-1 high, margins 7.5 bits; db).** A homothallic species with MAT1-1-1 and MAT1-2-1 fragments at 26-31% coverage. The new `high` MAT1-1 is not supported. It is caught by the margin rule.
3. **The #137 margin fallback, now measured campaign-wide: 380 `high` -> `medium` and 3,736 loci with a single-candidate identity margin.** *Macrophomina phaseolina* 228 drop because the margin moved from 5.18 to 4.96. This is a threshold effect from a db-driven 0.2-point identity shift, not evidence against MAT1-2. Decide the reporting rule before the site is rebuilt.
4. ***Pyricularia oryzae*** **93 `high` -> `medium` with a larger margin.** A distant (36 kb) unpolished MAT1-1-2 fragment enters the cluster, and the overlap loser MAT1-1-3 is no longer counted in `confidence_ignored_genes`. The tier drops and the span grows to 51 kb (db-driven). This is the #129 distant-fragment mechanism working in reverse.
5. **5 Hypocreales relaxed `idiomorph_gene_only` MAT1-1 calls lost**, db-driven in 5/5 re-runs: *Cordyceps* GCA_030783945 and YNP15, *Ophiocordyceps* YNP83, *Simplicillium* YNP88, *Clonostachys* GCA_041380665. The old calls had MAT1-1-1 at 41-62% identity (123-176 bits) and a polished MAT1-1-2. The new runs withhold every cluster by `modelled_gene_bar`. These look like real loci, lost to a record change. Also *Neopestalotiopsis* x3 (MAT1-2-1 at 31.5%, 72 bits, margin 33.9) became uncalled.
6. **Unlinked or second-locus homothallics lose the homothallic hint (11 genomes):** *Neosartorya* sp. G127-G150_14 x5, *A. hiratsukae* x2, *A. spinosus*, *A. latus* (an allodiploid hybrid), *A. felis* and 3 *A. fumigatus*. The MAT1-1 locus is now `mat_locus`, which is correct for the in-locus fragment, but each genome also has a full MAT1-2 at another locus. Genome-level reporting should show both.
7. ***Komagataella*** **(160 NEW + 15 locus count): MATalpha is never called.** It is withheld (`modelled_gene_bar`) in every genome. This is the known recall gap from the curation note: the GS115 MATalpha1 is `present: false`, and there is no CBS 7435 MATalpha record. Curator question 2 of that note is still open.
8. **One move from agree to disagree on the box-ORF scan:** Dothideomycetes sp. NU30, GCA_026013035.2, MAT1-2 high -> MAT1-1 high (db). The scan reads "MAT1-2: alpha-box ORF at the locus" (best alpha 186 bits, HMG 105). Both boxes are at the locus, so neither label is supported. *Setomelanomma holmii* and *Ophiobolus disseminans* flip the other way (MAT1-1 -> MAT1-2) inside the same scan class. All three are both-gene loci.
9. **Unverified homothallic conversions outside the literature-backed clades:** *Aspergillus alliaceus* x4 and 6 other *Aspergillus* species, *Didymella pinodes* 8, *Malbranchea*, *Onygena*, *Epibryon*. These rest on the policy's gene evidence only.
10. **Determinism, cosmetic:** the GAP1/BGL1 order inside `genes_found` differs between the pilot and the campaign in 21 *Pleosporaceae* genomes run with the same code apart from #140. The calls are identical; the list order is not stable.
11. ***Eutypa lata*** **12 LOST** (db). The old calls were 27%-identity relaxed MAT1-2-1 hits of 38-42 bits. One other *Eutypa* is a NEW call (needs new code + new db). These are weak either way, and there is no truth.

## What I could not verify

- There is no BioSample or read truth for *Parastagonospora*, *Macrophomina*, *Pyricularia*, *Aspergillus fumigatus* (beyond 2 genomes, unchanged), *Fusarium*, the yeasts, or the NEW Hypocreales/Xylariales calls. Those verdicts rest on report evidence, records and the Dothideomycetes box-ORF scan, which only covers Dothideomycetes.
- **Which db records cause each db-driven change.** Only the db-vs-code split was established. The *Macrophomina* identity shift, the *Pyricularia* MAT1-1-2 fragment and the Hypocreales losses are not traced to a specific record.
- The homothallic status of *A. alliaceus*, *D. pinodes*, *Stemphylium vesicarium*, *Phyllosticta capitalensis* and *Friedmanniomyces*. No literature was checked in this pass. The *F. graminearum* complex, Ptr, *S. sclerotiorum*, *Neurospora* FGSC and *Ceratocystis* rely on the pilot's literature checks.
- Sampling depth. The largest groups are covered by 3-9 re-runs each (main plus supplementary batch): *Macrophomina* 9 of 228, *Pyricularia* 9 of 191, *Fusarium* undetermined 8 of 430, *Parastagonospora* 8 of 87, *A. fumigatus* class 22 of 196. That 100% of each sample shares one cause, and every member of the group carries the same report signature (`attr.tsv`), supports extending the cause to the whole group. It is not a full census.
- No `nohomo` re-runs for the confidence and class categories. The single homothallic-linked DOWN call was explained by inspection.
- Per-locus `start`/`end` shifts within 1 kb were not judged: 2,925 genomes, median 6 bp.
