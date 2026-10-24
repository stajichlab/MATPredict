# The site-data step for the clean full-BFD run (2026-10-24)

Branch `site-data-cleanrun`. Nothing here deploys, tags or rebuilds the committed shards. Context: the clean-run plan
(`2026-10-22_clean-run-plan.md`, sections 8 and 10) and the curator's reporting rulings of 2026-10-22/24 (#129, #137,
#138; #145 is not covered).

## 1. What changed

| file | change |
|---|---|
| `scripts/build_cleanrun_tables.py` (new) | Per campaign: reads the per-genome reports and the v0.6.0 baseline reports, writes `genomes.tsv`, `loci.tsv`, `by_<rank>.tsv`. The columns written by v0.6.0 `analyze_asco.py` come first, unchanged. New columns: `margin_basis`, `baseline_confidence`, `tier_note`, `label_unsupported`, `high_basis` (per locus) and `genus`, `baseline_status`, `change_vs_baseline`, `n_tier_note`, `n_label_unsupported`, `high_basis` (per genome). The rules are in the module docstring. Reports are never changed. |
| `scripts/build_site_data.py` | The result folders are set in one place, `RELEASES` (`cleanrun` is the default; `v060` is still selectable). `--asco/--basidio/--mucoro` override one folder. A missing folder stops the build. Reads both table generations. For Mucoromycotina it reads the new `genomes.tsv`/`loci.tsv` (with `source`) or the v0.6.0 `calls.tsv`. Passes the new columns into the site data (below). |
| `scripts/build_locus_shards.py` | Uses the same `RELEASES`/overrides. Finds reports at any depth in `reports_all.tar.zst`. Appends `tier_note`, `label_unsupported`, `high_basis` to each locus when the campaign folder holds a clean-run `loci.tsv`. Duplicate run folders follow the table builder's rule. |
| `web/intro/landing.html` | The code label in the notice comes from `summary.json` (`release.code`). The old text is the fallback. |
| `web/intro/README.md` | Points here. |
| tests | `tests/test_build_cleanrun_tables_script.py` (23), `tests/test_build_site_data_script.py` (6). |

**How the new columns reach the site, without breaking current readers:**
- `genomes_<ds>.tsv` gets four appended columns: `label_unsupported`, `tier_notes`, `high_basis` and `locus_notes`. In `locus_notes`, each locus is `-` or its non-empty values joined by `,`, `;`-separated in `loci` order.
- `summary.json` gets per-dataset `label_unsupported`, `tier_note` and `high_basis` counts, plus a top-level `release` block.
- `taxa.json` gets a 15th row field, `unsupported`: genomes with at least one marked locus. `charts.js` reads only fields 0-13.
- The shards get fields 7-9. `locus-view.js` reads fields 0-6.

No page displays the new fields yet (see open questions).

**Checked:**
- A v0.6.0 build with the new script equals the old script's output once the appended fields are dropped (`taxa.json` rows, `summary.json` datasets, every genome TSV).
- The v0.6.0 shards rebuilt by the new script have the same content as the committed ones. Only the gzip header's OS byte differs, which comes from the Python version.
- A local preview (clean Ascomycota tables, v0.6.0 Basidiomycota and Mucoromycotina) loads with 0 console errors and 0 failed requests. The landing cards and charts render. `results.html` Aureobasidium drill-down draws 44 locus maps from the new shard.

## 2. Acceptance on the Ascomycota campaign (`2026-10-23_ascomycota_cleanrun` vs v0.6.0)

| check | judgement (ad-hoc scripts) | `build_cleanrun_tables.py` |
|---|---|---|
| legacy columns of `genomes.tsv` / `loci.tsv` | `analyze_asco.py` on HPCC | identical in all 19,415 / 20,985 rows |
| genome categories (`all_genomes_final.tsv`) | 14 categories | identical for all 19,415 genomes |
| 588 confidence-DOWN genomes | 380 #137 / 93 *Pyricularia* / 55 MTLA1 / 30 *Aureobasidium* | 380 `margin_fallback` / 100 `distant_fragment` / 55 `polished_model_lost` / 32 `vote_collapse` / 21 no note. The judgement's own mechanism table has the same numbers (A 380, D 100, C 55, B 32, class changed 16, other 5); the 93 and 30 are the *Pyricularia* and *Aureobasidium* parts of D and B. The rest of D is *Pseudogymnoascus* 3, *Penicillium* 2, *Ophiocordyceps* 1 and *Ambrosiella* 1; the rest of B is *Friedmanniomyces* 1 and *Phyllosticta* 1. |
| (same, per locus) | n/a | One genome, *Metschnikowia hawaiiensis* GCA_002374695.1, carries its MTLA1 drop on its second locus. The judgement's sorted-first pairing compared two different loci there; overlap matching puts the note on the locus that dropped. |
| label_unsupported | 757 (bit margin < 20 with >= 2 candidates, or *Aureobasidium*/*Xylaria*) | 781 = the same 757 loci + 24 *Friedmanniomyces* loci added by the curator's third genus (the other 25 of its 49 loci were already caught by the margin rule) |
| identity-fallback loci with margin < 20 | 3,736, not marked | none of the 3,736 is marked. 3,752 identity-fallback loci have a margin < 20; the extra 16 are *Aureobasidium*/*Xylaria* loci, which the genus rule marks (the judgement counted them under the genus rule too). |
| *Pleosporaceae* `high` split | 296 / 11 / 0 | 296 `supported` / 11 `tie_break` / 0 `grey`. The 5 `high` homothallic-candidate loci (`undetermined`) are blank: "the split does not apply". |
| campaign-wide (`--high-basis-scope all`) | 9,163 / 29 / 8 | 8,639 / 24 / 8 for idiomorph calls. The judgement also counted 529 `high` homothallic-candidate loci (524 supported, 5 tie-break): 8,639 + 524 = 9,163 and 24 + 5 = 29. |

**Tier drops campaign-wide**, all loci and not only the 588 primary loci:
- `margin_fallback` 381, `distant_fragment` 101, `polished_model_lost` 57, `vote_collapse` 51.
- 78 drops have no note (class changed, or no rule applies).

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
| per-dataset `code` | `scripts/build_site_data.py` (`RELEASES`) | same field |

**Still hand-written:**

| what | file:line | source to re-derive from |
|---|---|---|
| "Five datasets, one pipeline version" | `landing.html:51` | number of datasets in `summary.json` (5 unless LCG/Jena leave the site, see open questions) |
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
| dataset labels "(held out)" | `scripts/build_site_data.py` (`meta` in `main`), `web/intro/results.html:46` | only if LCG/Jena leave the site |

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

**2. Pack the reports for the shards** (per campaign). Run from inside the campaign folder:
```
find wave_* runs -name detection_report.yaml -o -name wall_seconds | sort \
  | tar -cf - -T - | zstd -19 -T8 -q -o reports_all.tar.zst
```
Use `wave_*` for Ascomycota, `runs` for Basidiomycota and Mucoromycotina.

**3. Commit the tables on a results branch.** Commit `README`, `MANIFEST.txt`, `lists/`, `genomes.tsv`, `loci.tsv`, `by_*.tsv` and `reports_all.tar.zst`; never the run folders. `.gitignore` ignores `results/*/runs/` only one level deep, so `wave_N/small_NNN/runs/` and `runs/small_NNN/runs/` are NOT ignored. Add the files by name, or extend `.gitignore` first.

**4. Set the release.** In `scripts/build_site_data.py` `RELEASES["cleanrun"]`:
- set the `mucoro` folder;
- set `code` to the release name;
- confirm the `asco`/`basidio` folders.

**5. Build the site data** (local, standard library only):
`python3 scripts/build_site_data.py`. This uses the default release, `cleanrun`, and writes `web/intro/data/{summary,taxa}.json` and the genome TSVs.

**6. Build the shards:** `python3 scripts/build_locus_shards.py` (needs PyYAML and `zstd`). Commit `web/intro/data/loci_*.json.gz`. These are the only generated files the Pages build does not regenerate.

**7. Edit the hand-written values** listed in section 3.

**8. Preview locally.** Run `python3 -m http.server --directory web/intro` and open `landing.html` and `results.html`. Check:
- the cards, the charts and a locus-map drill-down;
- the browser console.

The browser tool (`page/index.html`) exists only in the assembled site. To check the whole site, run `bash web/spike/page/build-site.sh <bundle> <site>`; it needs the WebAssembly bundle, which the CI builds.

**9. The curator decides; then tag.** A `v*` tag pushed to GitHub runs `.github/workflows/pages.yml`, which builds the bundle, runs `build-site.sh` (which re-runs `build_site_data.py` with the default release from the committed tables) and **deploys to the live Pages site, replacing it**. No archive of the old deployment was found. Keep the v0.6.0 data selectable (`--release v060`) for a rollback: re-tag with `RELEASES`' default set back.

**Verified:**
- The table builder reproduces the Ascomycota judgement (section 2).
- The v0.6.0 site data and shard contents are unchanged by the new code.
- The new and old Mucoromycotina table formats give the same site rows.
- The Basidiomycota v0.6.0 reports give the same site rows apart from 5 listed genomes with no report. The v0.6.0 Basidiomycota table left those 5 out; the new one lists them as `no_report`.
- The local preview loads with no errors.

**Not verified:**
- Whether a manual dispatch can publish to a preview location: **no.** `pages.yml` has `workflow_dispatch:` with no inputs and the same two jobs as a tag push. A dispatch builds from the chosen branch and deploys to the one `github-pages` environment, i.e. **it replaces the live site exactly like a tag**. A preview needs either a separate workflow (upload the artifact only, no `deploy` job) or a second Pages repository. Not run.
- **The Pages build fails once this branch merges, until the clean-run folders are committed.** The default release is `cleanrun`, and `build-site.sh` calls `build_site_data.py` without arguments. A failed build stops the deploy, so stale data is not published. To deploy before then, set `DEFAULT_RELEASE = "v060"`.
- The Basidiomycota clean run on the new layout: the layout is covered by a synthetic test (`runs/small_NNN/runs/<genome>`) and by the v0.6.0 tarball. It has not been run on the real folder, which is still being produced.
- The Mucoromycotina clean-run layout is assumed to repeat v0.6.0 (`runs/<SOURCE>__<id>`, `inputs.tsv`).
- The CI's Python 3.12 has not run the new `build_site_data.py`. It is standard library only and was run under the local `python3`.

## 5. Open questions for the curator

1. **LCG/Jena on the site.** Plan item 1b says the re-run held-out genomes "stay out of training and the site", but the current site shows both as held-out datasets. Keep them (the builder handles both), or drop them?
2. **Tier drops with no note (78 loci; 21 of the 588 genomes):**
   - 16 genomes where the locus class also changed, e.g. *Fusarium* `mat_locus` -> `partial_locus`;
   - 5 with no matching rule.
   The ruling gives four values and "blank when the tier did not drop", so these are blank with `baseline_confidence` showing the drop. A fifth value (`class_changed` / `other`)?
3. ***Friedmanniomyces* genus rule.** It adds 24 loci beyond the margin rule: 21 with a bit margin >= 20, 2 with no margin and 1 with an identity-fallback margin of 26.4 (judgement item 6.1: "add it, or rely on the margin rule"). Confirm.
4. **high_basis scope.** It is set by the BFD `FAMILY` = *Pleosporaceae* (as in the judgement), not by the scoped roster's taxid scope. The homothallic `high` loci are left blank. Confirm both.
5. **Showing the columns.** Should the explorer or landing page show `label_unsupported` / `tier_note` / `high_basis` (e.g. a mark on locus cards, a count per taxon), or is publishing them in the data and downloads enough for this release?
