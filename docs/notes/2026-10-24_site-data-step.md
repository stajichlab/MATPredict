# The site-data step for the clean full-BFD run (2026-10-24)

Branch `site-data-cleanrun`. Nothing here deploys, tags or rebuilds the committed shards. Context: the clean-run plan
(`2026-10-22_clean-run-plan.md`, sections 8 and 10) and the curator's reporting rulings of 2026-10-22/24 (#129, #137,
#138; #145 is not covered).

## 1. What changed

| file | change |
|---|---|
| `scripts/build_cleanrun_tables.py` (new) | Per campaign: reads the per-genome reports and the v0.6.0 baseline reports, writes `genomes.tsv`, `loci.tsv`, `by_<rank>.tsv`. The columns written by v0.6.0 `analyze_asco.py` come first, unchanged. New per-locus columns: `margin_basis`, `baseline_confidence`, `tier_note`, `label_unsupported` (`yes`/`no`), `label_unsupported_reason`, `high_basis`, `high_basis_outside_scope`. New per-genome columns: `genus`, `baseline_status`, `change_vs_baseline`, `n_tier_drop`, `n_tier_note`, `n_label_unsupported`, `high_basis`. The rules are in the module docstring. Reports are never changed. A locus with a null span gets no baseline match instead of aborting the table. |
| `scripts/build_site_data.py` | The result folders are set in one place, `RELEASES`. **The default stays `v060`** until the results PR commits the clean-run folders and flips `DEFAULT_RELEASE`; `cleanrun` is selectable with `--release`. `--asco/--basidio/--mucoro` override one folder. A missing folder stops the build. **Each dataset carries the code of its own folder** (`MANIFEST.txt` `tag:`/`commit:`, else the release that names the folder, else `unrecorded`). A build whose folders carry different codes is refused unless `--allow-mixed`, so v0.6.0 data can never be labelled with the clean run's code. Reads both table generations; for Mucoromycotina, the new `genomes.tsv`/`loci.tsv` (with `source`) or the v0.6.0 `calls.tsv`. |
| `scripts/build_locus_shards.py` | Uses the same `RELEASES`/overrides. Finds reports at any depth in `reports_all.tar.zst`. Picks the run folder of a duplicated genome with the table builder's own `choose_run`. Appends `baseline_confidence`, `tier_note`, `label_unsupported_reason`, `high_basis` to each locus when the campaign folder holds a clean-run `loci.tsv`. Stops when `loci.tsv` and the reports disagree on a genome's number of loci. |
| `web/intro/landing.html` | The code label comes from `summary.json` (`release.code`; the old text is the fallback). **The headline cards count BFD only**, and a fourth card gives the held-out ZyGoLife genomes as a separate number. The results text says the held-out sets are separate. The allele chart's group is labelled "ZyGoLife LCG + Jena (held out)". |
| `web/intro/results.html` | The ZyGoLife checkbox is labelled "ZyGoLife LCG and Jena (held-out sets, not counted with BFD)" and is **off by default**, so the default explorer view is BFD only. |
| `web/intro/README.md` | Points here. |
| tests | `tests/test_build_cleanrun_tables_script.py` (27), `tests/test_build_site_data_script.py` (12). |

**How the new columns reach the site, without breaking current readers:**
- `genomes_<ds>.tsv` gets five appended columns: `n_label_unsupported` (a count, `0` for any genome when the tables have the columns, blank only for v0.6.0 tables), `tier_notes`, `high_basis`, `loci_baseline_confidence` and `locus_notes`.
  - `loci_baseline_confidence` is the v0.6.0 tier of each locus, or `-` where nothing matched. A drop with no `tier_note` is visible here next to the reported tier (curator 2026-10-24).
  - In `locus_notes`, each locus is `-` or its non-empty values joined by `,`.
  - Both are `;`-separated in `loci` order.
- `summary.json` gets per-dataset `code` (from the folder), `held_out`, `n_label_unsupported`, `tier_drop`, `tier_drop_without_note`, and the `tier_note`, `high_basis` and `high_basis_outside_scope` counts. At the top level it gets `release` (name, code, folders, per-folder codes) and `headline` (BFD datasets only).
  - Folders inside the repository are recorded relative to it; others are recorded as absolute paths.
- `taxa.json` gets two row fields, 14 and 15. `unsupported` is the number of genomes with at least one marked locus. `high_basis_outside_scope` is `{basis: loci}`. `charts.js` reads only fields 0-13.
- The shards get fields 7-10. `locus-view.js` reads fields 0-6.

**Held-out sets (curator 2026-10-24).** The current data already split them: separate datasets `zygolife_lcg` and `zygolife_jena`, labelled "(held out)", source `ZyGoLife`; taxa rows with dataset `ZyGoLife` and cohort `LCG`/`JENA`. That split is kept, and `held_out: true` is added. Before this change, the landing cards summed all five datasets, so the 685 held-out genomes were inside "genomes searched". Now they are a separate card, and `summary.json` `headline` is BFD-only.

No page displays the new fields yet (see open questions).

**Checked:**
- A v0.6.0 build with the new script equals the old script's output once the appended fields are dropped (`taxa.json` rows, `summary.json` datasets, every genome TSV).
- The v0.6.0 shards rebuilt by the new script have the same content as the committed ones. Only the gzip header's OS byte differs, which comes from the Python version.
- A local preview (clean Ascomycota tables, v0.6.0 Basidiomycota and Mucoromycotina, built with `--allow-mixed`) loads with 0 console errors and 0 failed requests:
  - the cards read "22,974 BFD genomes searched" and "685 held-out ZyGoLife genomes (LCG, Jena), shown separately";
  - the code label reads "mixed: asco cleanrun-2026-10b (56e11cd); basidio v0.6.0 ...";
  - the explorer opens with only BFD checked (22,974 genomes);
  - the *Aureobasidium* drill-down draws 44 locus maps from the new shard.
- Without `--allow-mixed`, the same build is refused.

## 2. Acceptance on the Ascomycota campaign (`2026-10-23_ascomycota_cleanrun` vs v0.6.0)

| check | judgement (ad-hoc scripts) | `build_cleanrun_tables.py` |
|---|---|---|
| legacy columns of `genomes.tsv` / `loci.tsv` | `analyze_asco.py` on HPCC | identical in all 19,415 / 20,985 rows |
| genome categories (`all_genomes_final.tsv`) | 14 categories | identical for all 19,415 genomes |
| 588 confidence-DOWN genomes | 380 #137 / 93 *Pyricularia* / 55 MTLA1 / 30 *Aureobasidium* | 380 `margin_fallback` / 100 `distant_fragment` / 55 `polished_model_lost` / 32 `vote_collapse` / 21 no note. The judgement's own mechanism table has the same numbers (A 380, D 100, C 55, B 32, class changed 16, other 5); the 93 and 30 are the *Pyricularia* and *Aureobasidium* parts of D and B. The rest of D is *Pseudogymnoascus* 3, *Penicillium* 2, *Ophiocordyceps* 1 and *Ambrosiella* 1; the rest of B is *Friedmanniomyces* 1 and *Phyllosticta* 1. |
| (same, per locus) | n/a | One genome, *Metschnikowia hawaiiensis* GCA_002374695.1, carries its MTLA1 drop on its second locus. The judgement's sorted-first pairing compared two different loci there; overlap matching puts the note on the locus that dropped. |
| label_unsupported | 757 (bit margin < 20 with >= 2 candidates, or *Aureobasidium*/*Xylaria*) | **759 = 757 - 22 + 24.** Minus 22: loci reporting `undetermined` are excluded (curator 2026-10-24: a locus that publishes no idiomorph has no label to support). They are 20 `homothallic_candidate` (*Neurospora* 14, *Xanthoria* 2, *Botryotinia*, *Malbranchea*, *Pyrenophora*, *Rusavskia*), 1 `idiomorph_gene_only` *Naumovozyma*, and 1 *Aureobasidium* `mat_locus` (GCA_043230705.1, marked by both rules before). Plus 24: *Friedmanniomyces* loci added by the curator's third genus; the other 25 of its 49 loci were already caught by the margin rule. Multiset check against the judgement's list: judgement minus script = exactly those 22; script minus judgement = the 24. Before the exclusion it was 781 (= 757 + 24). |
| identity-fallback loci with margin < 20 | 3,736, not marked | none of the 3,736 is marked. 3,752 identity-fallback loci have a margin < 20; the extra 16 are *Aureobasidium*/*Xylaria* loci, which the genus rule marks (the judgement counted them under the genus rule too). |
| *Pleosporaceae* `high` split | 296 / 11 / 0 | 296 `supported` / 11 `tie_break` / 0 `grey`. The 5 `high` homothallic-candidate loci (`undetermined`) are blank: "the split does not apply". The scope stays the BFD FAMILY label (curator 2026-10-24). |
| outside *Pleosporaceae* (separate count, `high_basis_outside_scope`) | 13 tie-break, 8 grey | 8,343 supported / **13 tie_break / 8 grey**. These are in `by_class.tsv` (`outside_scope_high_*`), in `summary.json` and in `taxa.json`, never in `high_basis`. |
| campaign-wide (`--high-basis-scope all`) | 9,163 / 29 / 8 | 8,639 / 24 / 8 for idiomorph calls (= 296 + 8,343, 11 + 13, 0 + 8). The judgement also counted 529 `high` homothallic-candidate loci (524 supported, 5 tie-break): 8,639 + 524 = 9,163 and 24 + 5 = 29. |

**Tier drops campaign-wide**, all loci and not only the 588 primary loci:
- `margin_fallback` 381, `distant_fragment` 101, `polished_model_lost` 57, `vote_collapse` 51.
- 78 drops have no note (class changed, or no rule applies). They stay blank (curator 2026-10-24). A blank `tier_note` means "no listed mechanism", not "no drop"; the drop shows in `baseline_confidence` (loci), `n_tier_drop` (genomes), `loci_baseline_confidence` (site genome table) and `tier_drop_without_note` (`summary.json`: 78 for Ascomycota).
- `change_vs_baseline` follows the judgement's method: it compares the first locus of each sorted list. So 7 *F. graminearum* genomes (GCA_000966635.1, GCA_025427435.1, GCA_013618305.1, GCA_018345315.1, GCA_018346775.1, GCA_018219475.1, GCA_018219485.1) are `confidence DOWN` although no locus dropped. Their old primary MAT1-1 `high` became `undetermined`, and the medium `partial_locus` sorted to the front. Read "588" as "588 by the judgement's method". Per locus, `n_tier_note` and `n_tier_drop` are 0 for them.

**The margin_fallback rule meets the curator's "crossed the threshold" wording on all 381 loci.** On the v0.6.0 locus:
- the margin was >= 5 on 371 of them (270 identity fallback, 101 bit margin);
- there was no margin on the other 10.

So "the v0.6.0 margin was >= 5" was not added as a condition: it would change nothing on this campaign.

**Determinism:**
- Two runs give byte-identical output.
- Reading the reports from a tarball or from a folder gives the same bytes.
- `--jobs 8` and `--jobs 1` give the same bytes.
- Full campaign: about 40 s on a laptop.

## 3. Hand-written values that must change at release

Line numbers are for this branch.

**Data-driven now:**

| what | file:line | how derived |
|---|---|---|
| code label in the notice ("v0.6.0, October 2026") | `web/intro/landing.html:19` (fallback text) | `summary.json` `release.code`, from `RELEASES[...]["code"]` in `scripts/build_site_data.py`. **Set the `cleanrun` code to the release version (e.g. `v0.7.0`) before the tag.** |
| per-dataset `code` | each result folder's `MANIFEST.txt` | `tag:`/`commit:` lines; a v0.6.0 folder gets `RELEASES["v060"]["code"]` |
| headline cards (BFD genomes, call rate, loci; held-out count) | `landing.html` script | `summary.json` datasets with `source == "BFD"`, and the others |

**Still hand-written:**

| what | file:line | source to re-derive from |
|---|---|---|
| "Three BFD datasets and two held-out ZyGoLife sets, one pipeline version" | `landing.html:51` | the datasets in `summary.json` |
| "such as Orbiliomycetes and Dipodascomycetes" as low-rate clades | `landing.html:70` | the call-rate chart. The subphylum guard and the `MATdc`/`MATyl` families change Dipodascomycetes; re-read the chart. |
| "at most 6 loci per family" | `landing.html:29` | code constant (polish cap), unchanged unless the code changed |
| "about 100 bits", "40% or more" | `landing.html:40` | MAT-gene gate constants |
| "39 bits", "52 of 52", "7 of 8" | `landing.html:41` | flank-carried rule audit note |
| "25 bits" margin floor | `landing.html:44` | Mucoromycota classifier `min_margin` |
| Zygo 23: "23/23 locus, 23/23 idiomorph" | `landing.html:75` | clean-run regression gate: Zygo 23 23/23 (plan section 10) |
| Jena "64", "61 called", "63/63" | `landing.html:76` | clean-run Mucoromycotina tables (`summary.json` zygolife_jena) |
| LCG "621", "548 called (88%)" | `landing.html:77` | same, zygolife_lcg |
| classifier typing "540", "0 wrong", "96/108 vs 9/189" | `landing.html:78` | classifier analysis note. Unchanged unless re-measured. |
| "Regression panels (v0.6.0)", "166 + 23", "0 changed loci" | `landing.html:79` | **must change**: the clean run's gate was the 163-genome panel + Zygo 23, with the changes traced to commits (`results/2026-10_cleanrun_regression/`), not "0 changed". The curator words this. |
| *A. fumigatus* `detect` "297", "288 agree (97.0%)", "8 assemblies holding both" | `landing.html:80`, `:89` | **likely to change**: the clean run turned 196 *Aspergillus* `homothallic_candidate` into `mat_locus` MAT1-1. Re-run the 297-assembly comparison on the clean-run calls. |
| Fola "148", "147/148 (99.3%)", "3/3", "10%" | `landing.html:88` | reads note. Not affected by the run. |
| *A. fumigatus* reads "331", "293/296 (99.0%)", "Seven of eight" | `landing.html:89` | reads note. The `detect` part ("288/297") is as above. |
| Coccidioides "558" | `landing.html:91` | fallback only. JS overwrites it from `cocci_reads.json`. The `detect_call` column in `datasets/reads/coccidioides/results/cocci2_reads_vs_detect.tsv` holds **v0.6.0 calls** and needs re-joining to the clean-run calls if the agreement line is to describe the new release. |
| caveat numbers (50%, 100-base, one base in 40) | `landing.html:103` | reads notes. Unchanged. |
| "eight *A. fumigatus* strains" | `landing.html:106` | reads note. Unchanged. |
| "run with MATPredict v0.6.0 in October 2026" | `web/intro/about.md:51` | release version |
| software version "0.6.3" (citation, BibTeX) | `web/intro/about.md:28`, `:34` | the release version; also `CITATION.cff` and `pyproject.toml`. The reports say `0.6.1` whatever the code is (#139). |
| "about 1.6 MB for all 27.8k loci" | `web/intro/README.md:16` | the shard sizes printed by `build_locus_shards.py` |
| "about 1.9 MB for all 27.8k loci" | `scripts/build_locus_shards.py:10` (docstring) | same |
| dataset labels "(held out)" | `scripts/build_site_data.py` (`meta` in `main`), `web/intro/results.html:46`, `landing.html` (allele chart group) | fixed text; the held-out sets stay on the site as a separate category (curator 2026-10-24) |

Not hand-written: the landing cards, the call-rate, confidence, locus-class and allele charts, the explorer, the genome tables and the locus maps. All of them come from `summary.json`, `taxa.json` and the shards.

## 4. Plan: the ordered site update (for after phase 5 of every campaign)

Paths are on the HPCC unless marked local. Steps 1-3 run once per campaign. Nothing is committed to `main` or tagged without the curator.

**1. Build the tables** (repeat per campaign; output into the campaign folder, next to `MANIFEST.txt`):
```
# Ascomycota
python3 scripts/build_cleanrun_tables.py --campaign results/2026-10-23_ascomycota_cleanrun \
  --baseline results/2026-10-03_ascomycota_v060/reports_all.tar.zst \
  --samples /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv \
  --out results/2026-10-23_ascomycota_cleanrun --jobs 8
# Basidiomycota (after jobs 30336695-30336710 finish and phase 5 is judged)
python3 scripts/build_cleanrun_tables.py --campaign results/2026-10-23_basidiomycota_cleanrun \
  --baseline results/2026-10-03_basidiomycota_v060/reports_all.tar.zst \
  --samples .../samples.csv --rank order --out results/2026-10-23_basidiomycota_cleanrun --jobs 8
# Mucoromycotina (folder name not yet known)
python3 scripts/build_cleanrun_tables.py --layout mucoro --campaign results/<mucoro_cleanrun> \
  --baseline results/2026-10-03_mucoromycotina_mat/reports_all.tar.zst \
  --baseline-calls results/2026-10-03_mucoromycotina_mat/calls.tsv --rank order --out results/<mucoro_cleanrun>
```
The Ascomycota `genomes.tsv`, `loci.tsv` and `by_class.tsv` replace the files `analyze_asco.py` wrote there. The old columns are byte-identical, so nothing is lost.

**2. Pack the reports for the shards** (per campaign). Run from inside the campaign folder, naming the one tree that exists:
```
# Ascomycota
find wave_* \( -name detection_report.yaml -o -name wall_seconds \) | sort | tar -cf - -T - | zstd -19 -T8 -q -o reports_all.tar.zst
# Basidiomycota and Mucoromycotina
find runs \( -name detection_report.yaml -o -name wall_seconds \) | sort | tar -cf - -T - | zstd -19 -T8 -q -o reports_all.tar.zst
```

**3. Commit the tables on a results branch.** Commit only the top-level files of each campaign folder: `README`, `MANIFEST.txt`, `lists/`, `genomes.tsv`, `loci.tsv`, `by_*.tsv` and `reports_all.tar.zst`. Never commit run folders. How `.gitignore` treats them (checked with `git check-ignore`):
- `results/*/runs/` (`.gitignore:13`) ignores the whole Basidiomycota and Mucoromycotina `runs/` tree.
- The Ascomycota `wave_N/small_NNN/runs/` is NOT ignored, because the pattern matches a `runs/` directly under the campaign folder only. Do not `git add` the campaign folder as a whole; add the files above by name (or add `results/*/wave_*/` to `.gitignore` in the results PR).
- The tables and `reports_all.tar.zst` sit at the top of the folder, outside both trees, so they are always added normally.

**4. Set the release** (in the results PR). In `scripts/build_site_data.py`:
- in `RELEASES["cleanrun"]`, set the `mucoro` folder, set `code` to the release name, and confirm the `asco`/`basidio` folders;
- set `DEFAULT_RELEASE = "cleanrun"`.

Every folder needs a `MANIFEST.txt` with the same `tag:` line. Otherwise the build refuses the mix. That is wanted: do not pass `--allow-mixed` for a release.

**5. Build the site data** (local, standard library only):
`python3 scripts/build_site_data.py`. With no arguments it exercises exactly what the Pages build will run. It writes `web/intro/data/{summary,taxa}.json` and the genome TSVs.

**6. Build the shards:** `python3 scripts/build_locus_shards.py` (needs PyYAML and `zstd`). Commit `web/intro/data/loci_*.json.gz`. These are the only generated files the Pages build does not regenerate.

**7. Edit the hand-written values** listed in section 3.

**8. Preview locally.** Run `python3 -m http.server --directory web/intro` and open `landing.html` and `results.html`. Check:
- the cards, the charts and a locus-map drill-down;
- the browser console.

The browser tool (`page/index.html`) exists only in the assembled site. To check the whole site, run `bash web/spike/page/build-site.sh <bundle> <site>`; it needs the WebAssembly bundle, which the CI builds.

**9. The curator decides; then tag.** A `v*` tag pushed to GitHub runs `.github/workflows/pages.yml`, which:
- builds the bundle;
- runs `build-site.sh`, which re-runs `build_site_data.py` with the default release from the committed tables;
- **deploys to the live Pages site, replacing it**.

No archive of the old deployment was found. For a rollback, set `DEFAULT_RELEASE` back to `v060`, commit, and re-tag. The v0.6.0 folders stay in `RELEASES`. The v0.6.0 shards are in git history; rebuild them with `build_locus_shards.py --release v060`.

**Verified:**
- The table builder reproduces the Ascomycota judgement (section 2).
- The v0.6.0 site data and shard contents are unchanged by the new code.
- The new and old Mucoromycotina table formats give the same site rows.
- The Basidiomycota v0.6.0 reports give the same site rows apart from 5 listed genomes with no report. The v0.6.0 Basidiomycota table left those 5 out; the new one lists them as `no_report`.
- The local preview loads with no errors.

**Not verified:**
- Whether a manual dispatch can publish to a preview location: **no.** `pages.yml` has `workflow_dispatch:` with no inputs and the same two jobs as a tag push. A dispatch builds from the chosen branch and deploys to the one `github-pages` environment, i.e. **it replaces the live site exactly like a tag**. A preview needs either a separate workflow (upload the artifact only, no `deploy` job) or a second Pages repository. Not run.
- **Merging this branch does not change the live data.** `DEFAULT_RELEASE` stays `v060`, so a `v*` tag before the results PR rebuilds the v0.6.0 site, now with the BFD-only headline cards and the explorer's ZyGoLife box unchecked.
- The Basidiomycota clean run on the new layout: the layout is covered by a synthetic test (`runs/small_NNN/runs/<genome>`) and by the v0.6.0 tarball. It has not been run on the real folder, which is still being produced.
- The Mucoromycotina clean-run layout is assumed to repeat v0.6.0 (`runs/<SOURCE>__<id>`, `inputs.tsv`).
- The CI's Python 3.12 has not run the new `build_site_data.py`. It is standard library only and was run under the local `python3`.

## 5. Curator rulings (2026-10-24) and what is still open

**Ruled and implemented:**
1. **LCG/Jena** stay on the site as a separate, clearly labelled held-out category, out of the BFD headline counts (section 1).
2. **Tier drops with no listed mechanism** (78 loci; 21 of the 588 genomes: 16 class also changed, 5 no rule) stay blank. `baseline_confidence` is shown next to the reported tier in the loci table, the genome tables (`n_tier_drop`, site `loci_baseline_confidence`), `summary.json` (`tier_drop`, `tier_drop_without_note`) and the shards.
3. ***Friedmanniomyces*** is on the genus list. It adds 24 loci beyond the margin rule: 21 with a bit margin >= 20, 2 with no margin, and 1 with an identity-fallback margin of 26.4.
4. **high_basis scope** stays the BFD FAMILY label *Pleosporaceae*, and homothallic `high` loci are blank. `high` calls outside the scope are counted separately (`high_basis_outside_scope`: 8,343 / 13 / 8 in Ascomycota).
5. **`undetermined` loci are never marked** `label_unsupported`, by either rule. This removes 22 loci: 781 -> 759.

**Still open:**
1. **Showing the columns.** Should the explorer or landing page show `label_unsupported`, `tier_note` or `high_basis` (e.g. a mark on locus cards, a count per taxon)? Or is publishing them in the data and downloads enough for this release?
2. **The 72 suppressed genomes** are in the clean run's `lists/` and appear as `no_report` (44 of them had v0.6.0 calls), counted in "BFD genomes searched". This is the same behaviour as v0.6.0. Should they be dropped from the lists, or given a `suppressed` status the site can hide?
3. **The hand-written values** in section 3 that need wording rather than a script: the regression-panel row, the *A. fumigatus* `detect` row, the Coccidioides `detect_call` join, the `about.md` versions, and the release code label.
