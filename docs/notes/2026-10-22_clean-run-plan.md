# Plan: the clean full-BFD run (draft, 2026-10-22)

Status: **draft for the curator; nothing has been run.** It follows the decisions recorded in `docs/notes/2026-10-19_session-handoff.md` (section 5) and the
curator's rulings of 2026-10-22. Facts come from the v0.6.0 campaign records (`results/2026-10-03_*_v060/README.md`, `docs/regression-check.md`) and the repository at `origin/main`
2d14712 plus the open PRs named below. Items I could not verify are marked **(unverified)**; open decisions are in section 9.

## 1. Purpose and scope
Re-run detection on the BFD genomes with the current code and database and replace the site data, **after** the new results have been diffed against v0.6.0 and the curator has reviewed the differences.
The old site data stay in place until then (curator ruling).
v0.6.0 baseline (three campaigns, all from frozen worktree `run-7c7ed99`):

| campaign | genomes | how | wall |
|---|---|---|---|
| Ascomycota | 19,415 | 9 waves of about 2,158, 64-way, partition `exfab`, `GENOME_TIMEOUT 3600`; pilot of 246 first (mean 158 s/genome, max 624 s) | about 1.5 h per wave of compute, plus queueing |
| Basidiomycota | 3,270 rows in `genomes.tsv`; the job lists total 3,275 (3 small jobs of 1,073 + 1,072 + 1,072 genomes under 500 Mb, 64-way, timeout 3600; 1 large job of 58 genomes at or above 500 Mb, 32-way, timeout 14400, account `exfab`); **the difference of 5 is unexplained (unverified)** | same worker pattern | |
| Mucoromycotina | 974 genomes in `inputs.tsv`: **BFD Mucorales + Umbelopsidales 289 (288 after the suppress list)**, LCG 621 and Jena 64 (the last two are held-out sets, never used for training, curation or classifier builds); `detect --phylum Mucoromycota --emit-cds-fasta`, **no taxid** (LCG and Jena have none), 4-chunk Slurm array 29371553 (`run_chunk.slurm`, `chunk_{0..3}.tsv`, committed), tables from `collect.py`; 47.4 CPU-h, median 183 s per genome | |
`scripts/build_site_data.py` reads all three (the *Mucoromycotina* part from `calls.tsv`, BFD rows only), so replacing the site's data needs all three. **Decided (curator, 2026-10-22): all three campaigns.** The *Mucoromycotina* production is documented in `results/2026-10-03_mucoromycotina_mat/NOTE.md` and its scripts are committed, so it is a small, reproducible add-on (288 BFD genomes, about 15 CPU-hours); whether to also re-run the 685 held-out LCG/Jena genomes is a separate question (section 9, decision 1b).

## 2. Gates (curator ruling 2026-10-22)
The run starts when these are merged:
1. #119 (reference exon structures): **merged**. #123 carries the three `record_version` bumps it missed.
2. #124 (`homothallic_candidate` reports idiomorph `undetermined`): **merged**.
3. The *Cochliobolus* scoped roster (**decided 2026-10-22: it joins, `--scoped-roster-policy on`**): the mechanism (a per-gene `taxonomic_scope`, flag default off), GAP1/BGL1 as optional one-member slots with a 50% identity floor, *Pleosporaceae* only. It is gated on its PR being merged after its Fable-review fixes (version bump of the C4/C5 reference records, loud failure of the scope lookup, corrected tier wording, documented BGL1 floor margin) and on the curator's re-review of the `5016_c4` / `5016_c5` records. Measured: 90 of 101 uncalled clade genomes recovered; 2,252 of 2,252 Dothideomycetes reports identical outside the scope. The floor must be re-checked before the scope is widened to any other family.
**Not gating:** PR A (more exon structures plus a validator rule), PR B (the *liq* / *C. auris* / *Okayama* curation corrections), roster growth (a later release if it proceeds), and the gene-model extension (shadow mode changes no call). PR B changes two *Teloschistales* and one *Serinales* reference protein, so it rides a later release.

## 3. Run configuration
- **Code:** a frozen worktree at a commit of `main` pinned by a non-`v` tag such as `cleanrun-2026-10` (never edited while jobs run), recorded in the manifest (section 7). `pyproject.toml` says 0.6.3 (last `v` tag v0.6.2); the version is not bumped for the run (**decision 3**).
- **Settings, passed explicitly (code defaults are unchanged):** `--partial-core-policy coverage_profile`, `--confidence-cap-policy off`, `--strong-single-gene-policy off`, `--homothallic-policy both_genes`, `--scoped-roster-policy on`; `--polish-core-pair` is the default (on); lineage routing from the taxid (default); no `--helixerlite`, no tracks.
  The companion-gene veto is active with the partial-core rule. The policy can act only on the Ascomycota MAT family **(verified in the code, 2026-10-22)**: `partial_core.load_core_profiles` returns `None` unless `db/<phylum>/classifiers/MAT_core/manifest.yaml` exists, which only Ascomycota has, and the homothallic evaluation returns the call unchanged when the profiles are absent; Basidiomycota (`homothallic_screen: false` at phylum level) and *Mucoromycotina* are therefore unaffected by construction.
- **Genetic code 26 defect of v0.6.0 is fixed** (exonerate longhand code; 20 genomes had failed).
- **Timeouts:** `GENOME_TIMEOUT` 3600 s, 14400 s for genomes at or above 500 Mb.
- **Tooling: committed and sufficient (corrected 2026-10-22; an earlier draft of this note wrongly called it a gap).** The v0.6.0 campaigns ran through `scripts/submit_binned_panel.py` (splits a sample list by genome size, under or at/above 500 Mb, and submits one `scripts/run_clade_panel.slurm` job per chunk with its own `OUT_DIR`, `rollout_summary.yaml` and `reports.tar.zst`; large genomes use `GENOME_TIMEOUT=14400`) and `scripts/run_clade_panel.slurm`, which runs `matpredict detect --genome ... --taxid ...` per genome (resumable: a genome with a non-empty `detection_report.yaml` is skipped; taxid routing and genetic code resolved per genome, warmed up once serially) and passes `${DETECT_ARGS}` straight to `detect`. `submit_binned_panel.py` submits with `--export=ALL,...`, so `DETECT_ARGS` set in the submitting shell reaches every job; the job log prints `detect_args=...`. The batch driver `scripts/run_detection_batch.py` / `batch_runner.run_batch` is a **different** path that does not resolve the genetic code per genome and forwards none of the newer policy flags; **it is not used for the campaign.**
  The run is therefore `DETECT_ARGS="--partial-core-policy coverage_profile --confidence-cap-policy off --strong-single-gene-policy off --homothallic-policy both_genes --scoped-roster-policy on" python3 scripts/submit_binned_panel.py --sample-list LIST.tsv --out <results dir> --worktree <frozen worktree> --name <name>` (run `--dry-run` first and check that the printed `sbatch` command and the job log carry the flags).
  **Provenance:** every `detection_report.yaml` records its parameters under `run.parameters` (`partial_core_policy`, `homothallic_policy`, `confidence_cap_policy`, `strong_single_gene_policy`, `polish_core_pair`, the polish cap and identity), so the settings of each genome can be read back from its report; the manifest of section 7 adds the code commit, database hash and environment hashes.
  **Still to check in the dry run:** that the Ascomycota 9-wave layout of v0.6.0 (`jobs.tsv` names `wave_0..8`, about 2,158 genomes each, classes and orders interleaved) can be reproduced with `submit_binned_panel.py`'s size bins (it names chunks `small_000...`), or that a sample list pre-split into interleaved waves is submitted one wave at a time; and the *Mucoromycotina* path (`results/2026-10-03_mucoromycotina_mat/run_chunk.slurm`: `--phylum Mucoromycota --emit-cds-fasta`, no taxid), which is a separate script.

## 4. Phases
| phase | what | done when |
|---|---|---|
| 0 | prerequisites: gates in section 2; frozen worktree and the non-`v` tag; `pixi.lock` hash; database content hash; `submit_binned_panel.py --dry-run` with `DETECT_ARGS` set; a 3-genome live test of the worker | the manifest exists, the printed `sbatch` command and the job log carry the flags, and the 3-genome test reproduces a known call |
| 1 | **regression gate** (`docs/regression-check.md`): `testset/regression_panel.tsv` plus Zygo 23, baseline versus candidate, partition `exfab`; every changed call listed and judged | no unexplained change; every change in an expected category (section 6) |
| 2 | **determinism check:** about 200 genomes run twice with the same code; diff | the noise level is known before differences are attributed to code (the handoff recorded one *Ogataea* `MATsc` start of 199 versus 349 between runs, uninvestigated) |
| 3 | **pilot** (list built: `results/2026-10-22_clean_run_pilot/`, 569 Ascomycota/Basidiomycota genomes plus 25 for the *Mucoromycotina* path, seed 20261022) of about 500 genomes: stratified by class like the 246-genome v0.6.0 pilot, plus the known problem clades (all 101 uncalled *Cochliobolus*/*Bipolaris*, *Pyrenophora*, *Aureobasidium*, *Botrytis*, *Hypocreales* with the *Beauveria*-type structure) | per-class runtime and timeout rate known; no systematic failure; the homothallic policy behaves as measured |
| 4 | **waves:** Ascomycota in 9 waves of about 2,158, 64-way on `exfab`; Basidiomycota small (3 jobs) and large (1 job, 32-way); *Mucoromycotina* (one 4-chunk array, `run_chunk.slurm` pattern, no taxid); failures re-run in a separate folder | at least 99% of genomes have a report; failures listed with causes |
| 5 | **diff against v0.6.0** (section 6) and curator review | every changed call explained or put on the curator's list |
| 6 | **site** (section 8), only on the curator's decision | curator approves; tag |

## 5. Resources and timing (estimates, from v0.6.0)
- Compute: 19,415 genomes x 158 s (v0.6.0 pilot mean) is about 3.1 million core-seconds, roughly 850 core-hours for Ascomycota; at 64-way a wave of 2,158 is about 1.5 h of compute. The new rules add per-genome work (partial-core, the homothallic evaluation, companion veto) so **the pilot must re-measure**; take the pilot mean times 1.3 as the planning figure until then.
- `exfab` is one 256-core node (checked 2026-10-22: 18 cores allocated, 238 idle, 3 queued jobs); it shares an account with unrelated jobs (`nf-*`, `nxf_*`), so queue time is unpredictable.
- Output volume: v0.6.0 kept a `reports_all.tar.zst` and tables per campaign; the same layout is proposed (`results/<date>_<phylum>_<tag>/`), with per-genome reports tarred by wave. Disk quota on `/bigdata` is **unverified**; check before wave 0.

## 6. The diff and its pre-registered expectations
Compare per genome and per locus against `results/2026-10-03_{ascomycota,basidiomycota}_v060/{genomes,loci}.tsv`: status (called or not), number of loci, families called, idiomorph, confidence, `locus_class`, gene models (spans), and wall time.
**Expected change categories** (merged work since v0.6.0, from the 2026-10-19 handoff; each must be counted, not assumed):
1. partial-core rule on (#95) with the companion-gene veto (#102): flips towards the profile call in Sclerotiniaceae/Botrytis-type loci; the flip review found 41 of 46 flips right (dominated by *Botrytis*);
2. polish-core-pair fix (#104) and `Ascomycota` core profiles (#86): gene-model and vote changes where a sexM/sexP-style pair decided;
3. records: *Ceratocystis* (`unisexual_switching`), *Botrytis*, *Parastagonospora*, *Komagataella* (`MATkp`), *Ogataea* complex (`MATog`), Dipodascales (`MATyl`, `MATdc`), *Hemileia* `rustHD`, remnant status for seven records;
4. the Dipodascales fallback guard: **losses are expected**, at least 52 *Lipomycetales* and 20 *Dipodascales* genomes lose calls until they get their own families;
5. `--homothallic-policy both_genes`: about 99 loci in the 3,181-genome test set become `homothallic_candidate` (idiomorph `undetermined` once #124 merges), mostly *Dothideomycetes*; no control changes;
6. *Cochliobolus*/*Bipolaris*: recovery from the scoped roster (90 of 101 in the audit list); `high` calls are reported split into *supported* (a gene beyond the core gene and the scoped flanks is polished, or reaches 50% identity with 30% coverage) and *tie-break* (every such extra gene is unpolished and under 30% coverage, or is the opposite idiomorph's core gene; issue #129). The thresholds are checked on the pilot list before the full run;
7. reference exon structures (#119, PR A): **no detection change expected**, which is itself a check.
**Acceptance (to be fixed before the run):** (a) at least 99% of genomes yield a report; (b) the regression panel shows no change outside the expected categories; (c) every genome whose call changed is assigned to one expected category or listed for the curator, with counts per category; (d) the overall call rate does not fall except by the Dipodascales guard; (e) the determinism check (phase 2) bounds the unexplained noise, and any difference below that level is not attributed to code.

## 7. Reproducibility record
A `run_manifest.yaml` per campaign, written at submission and committed with the results: code commit and tag, worktree path, `pixi.lock` and environment hashes, database content hash (the registry hash of `db/`), exact flags, genome list and its source (`samples.csv` / `ncbi_accessions_taxonomy.csv` hashes), Slurm job ids, the node counts, and the dates. v0.6.0 recorded the frozen worktree but not the flags in one place; the manifest closes that.

## 8. The site
Per the 2026-10-19 handoff **(reconstructed, not run)**: `.github/workflows/pages.yml` deploys on `v*` tags and on manual dispatch; `scripts/build_site_data.py` reads **hard-coded** result folders (`2026-10-03_ascomycota_v060`, `..._basidiomycota_v060`, `..._mucoromycotina_mat`, and the *Coccidioides* reads table); locus-map shards `web/intro/data/loci_*.json.gz` are committed and rebuilt by `scripts/build_locus_shards.py`, not by the Pages build; some numbers in `web/intro/landing.html` are hand-written; **a tag overwrites the live site and no archive of old deployments was found**.
So replacing the data means: new dated result folders, edited paths in `build_site_data.py`, rebuilt shards, edited landing numbers, a local build, and only then a tag, on the curator's explicit decision. Whether a manual dispatch can publish to a preview location is **unverified**; check before relying on it.

## 9. Decisions for the curator

**All six are now decided (2026-10-22):** (1) all three campaigns, and (1b) the 685 held-out LCG/Jena genomes included; (2) order Ascomycota, Basidiomycota, *Mucoromycotina*; (3) a non-`v` run tag; (4) the determinism check; (5) no decision needed (Ascomycota-only by construction); (6) run without the *Cochliobolus* roster if it is not ready. Their text is kept below. **Item 2's order was recommended and not separately confirmed; treat it as the working order unless the curator changes it.**

1. **Scope: decided (2026-10-22): all three campaigns.** **1b: decided (2026-10-22): also re-run the 685 held-out LCG and Jena genomes** (about 35 CPU-hours at the v0.6.0 median of 183 s) with the same code, as the held-out set for roster and classifier work; they stay out of training and the site.
2. **Order:** run Ascomycota first (largest, most affected), then Basidiomycota, then *Mucoromycotina*? Recommendation: yes.
3. **Version and tag: decided (2026-10-22): a tag that does not start with `v`** (for example `cleanrun-2026-10`), so the run's commit is pinned by name without matching the Pages workflow's `v*` trigger (verified: `.github/workflows/pages.yml` runs on `tags: ["v*"]` and on manual dispatch, so any `v` tag, including `v0.7.0-rc1`, would deploy). A real `v` tag is created only when the curator decides to replace the site.
4. **Determinism check (phase 2): decided (2026-10-22): yes**, about 200 genomes run twice with identical code, before the pilot; the diff is the noise floor for section 6, and a cause for the *Ogataea* `MATsc` coordinate flip is looked for.
5. **Homothallic policy outside Dothideomycetes: resolved by the code, no decision needed.** It acts only on the Ascomycota MAT family (section 3); the pilot reports its count per class (it was measured on Dothideomycetes, 3,181 genomes), and no other phylum is affected.
6. **The *Cochliobolus* scoped roster is on in the run (decided 2026-10-22, superseding the earlier run-without-it fallback).** If it is not merged and re-reviewed when the run is ready, the fallback stands: run without it and, if it validates afterwards, re-run only the *Pleosporaceae* genomes as a separate dated folder and swap those rows in (the strong-single-gene check after the run is the other route to the same gain).

## 10. Risks
Queue contention (shared account); timeouts on large genomes (the 14400 s class); empty or missing FASTA (for example `GCA_986281625.1_DAH1041FM` has an empty file) listed rather than retried; a database or code change during the run (the frozen worktree prevents it); the intermittent test (`tests/detect/test_pipeline.py::test_fast_path_missing_gene_is_polished_in_the_existing_clusters_window`) also fails on unmodified main and is unrelated to the run; and the unverified items above.

## 11. After the run
Re-derive the family counts in `docs/notes/2026-10-18_families-needing-curation.md` from the new results (they come from v0.6.0); run the 12-locus sequence check for the strong-single-gene rule (decided to follow the clean run); re-run the Basidiomycota scope-only family table; update the failure tracker's revisit list with the new numbers.
