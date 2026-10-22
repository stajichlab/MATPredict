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
| Mucoromycotina | produced by a different table (`results/2026-10-03_mucoromycotina_mat/calls.tsv`; how it was produced is **unverified**) | | |
`scripts/build_site_data.py` reads all three, so a replacement of the site needs all three (or an explicit decision to keep one at v0.6.0).

## 2. Gates (curator ruling 2026-10-22)
The run starts when these are merged:
1. #119 (reference exon structures): **merged**. #123 carries the three `record_version` bumps it missed.
2. #124 (`homothallic_candidate` reports idiomorph `undetermined`): open.
3. The *Cochliobolus* decision: a scoped *Pleosporaceae* GAP1/BGL1 roster is being prototyped (branch `cochliobolus-scoped-roster`); its recovery and control results decide whether it joins.
**Not gating:** PR A (more exon structures plus a validator rule), PR B (the *liq* / *C. auris* / *Okayama* curation corrections), roster growth (a later release if it proceeds), and the gene-model extension (shadow mode changes no call). PR B changes two *Teloschistales* and one *Serinales* reference protein, so it rides a later release.

## 3. Run configuration
- **Code:** a frozen worktree at a tagged commit of `main` (never edited while jobs run), recorded in the manifest (section 7). Version: `pyproject.toml` says 0.6.3, last tag v0.6.2 (**decision 3**).
- **Settings, passed explicitly (code defaults are unchanged):** `--partial-core-policy coverage_profile`, `--confidence-cap-policy off`, `--strong-single-gene-policy off`, `--homothallic-policy both_genes`; `--polish-core-pair` is the default (on); lineage routing from the taxid (default); no `--helixerlite`, no tracks.
  The companion-gene veto is active with the partial-core rule. `homothallic_screen` applies only to families that declare it, so the policy affects Ascomycota MAT (**unverified for other phyla**: the pilot will show).
- **Genetic code 26 defect of v0.6.0 is fixed** (exonerate longhand code; 20 genomes had failed).
- **Timeouts:** `GENOME_TIMEOUT` 3600 s, 14400 s for genomes at or above 500 Mb.
- **Script gap (must be closed first).** The batch driver (`scripts/run_detection_batch.py`, `src/MATPredict/detect/batch_runner.py`) does **not** forward any of the newer policy flags (no `partial_core_policy`, `confidence_cap_policy`, `strong_single_gene_policy` or `homothallic_policy` in either file), and the v0.6.0 campaign worker is **not committed** in the repository. The 2026-10-18 homothallic runs used a per-genome `python -m MATPredict detect ...` worker (`results/2026-10-18_homothallic_patch/patched_worker.slurm`: resumable, per-genome timeout, fasta decompressed to `$SCRATCH`, report kept, a check that the code imported is the worktree's). The plan is to **commit a campaign worker** (`scripts/run_campaign_wave.slurm` plus a list builder) based on that pattern, with the flags in one place and a manifest written at submission. This is a small code PR and is a precondition, not an option.

## 4. Phases
| phase | what | done when |
|---|---|---|
| 0 | prerequisites: gates in section 2; the campaign-worker PR merged; frozen worktree and tag; `pixi.lock` hash; database content hash; a dry run of the worker on 3 genomes | the manifest exists and a 3-genome dry run reproduces a known call |
| 1 | **regression gate** (`docs/regression-check.md`): `testset/regression_panel.tsv` plus Zygo 23, baseline versus candidate, partition `exfab`; every changed call listed and judged | no unexplained change; every change in an expected category (section 6) |
| 2 | **determinism check:** about 200 genomes run twice with the same code; diff | the noise level is known before differences are attributed to code (the handoff recorded one *Ogataea* `MATsc` start of 199 versus 349 between runs, uninvestigated) |
| 3 | **pilot** of about 500 genomes: stratified by class like the 246-genome v0.6.0 pilot, plus the known problem clades (all 101 uncalled *Cochliobolus*/*Bipolaris*, *Pyrenophora*, *Aureobasidium*, *Botrytis*, *Hypocreales* with the *Beauveria*-type structure) | per-class runtime and timeout rate known; no systematic failure; the homothallic policy behaves as measured |
| 4 | **waves:** Ascomycota in 9 waves of about 2,158, 64-way on `exfab`; Basidiomycota small (3 jobs) and large (1 job, 32-way); Mucoromycotina as decided; failures re-run in a separate folder | at least 99% of genomes have a report; failures listed with causes |
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
6. *Cochliobolus*/*Bipolaris*: recovery if the scoped roster joins;
7. reference exon structures (#119, PR A): **no detection change expected**, which is itself a check.
**Acceptance (to be fixed before the run):** (a) at least 99% of genomes yield a report; (b) the regression panel shows no change outside the expected categories; (c) every genome whose call changed is assigned to one expected category or listed for the curator, with counts per category; (d) the overall call rate does not fall except by the Dipodascales guard; (e) the determinism check (phase 2) bounds the unexplained noise, and any difference below that level is not attributed to code.

## 7. Reproducibility record
A `run_manifest.yaml` per campaign, written at submission and committed with the results: code commit and tag, worktree path, `pixi.lock` and environment hashes, database content hash (the registry hash of `db/`), exact flags, genome list and its source (`samples.csv` / `ncbi_accessions_taxonomy.csv` hashes), Slurm job ids, the node counts, and the dates. v0.6.0 recorded the frozen worktree but not the flags in one place; the manifest closes that.

## 8. The site
Per the 2026-10-19 handoff **(reconstructed, not run)**: `.github/workflows/pages.yml` deploys on `v*` tags and on manual dispatch; `scripts/build_site_data.py` reads **hard-coded** result folders (`2026-10-03_ascomycota_v060`, `..._basidiomycota_v060`, `..._mucoromycotina_mat`, and the *Coccidioides* reads table); locus-map shards `web/intro/data/loci_*.json.gz` are committed and rebuilt by `scripts/build_locus_shards.py`, not by the Pages build; some numbers in `web/intro/landing.html` are hand-written; **a tag overwrites the live site and no archive of old deployments was found**.
So replacing the data means: new dated result folders, edited paths in `build_site_data.py`, rebuilt shards, edited landing numbers, a local build, and only then a tag, on the curator's explicit decision. Whether a manual dispatch can publish to a preview location is **unverified**; check before relying on it.

## 9. Decisions for the curator
1. **Scope:** Ascomycota only (the 19,415 of the handoff), or all three campaigns so the site can be replaced as a whole? Recommendation: Ascomycota and Basidiomycota; keep Mucoromycotina at v0.6.0 until its production is understood.
2. **Order:** run Ascomycota first (largest, most affected), then Basidiomycota? Recommendation: yes.
3. **Version and tag:** tag the run commit (for example `v0.7.0-rc1`) or run from an untagged frozen commit? Recommendation: a release-candidate tag (a tag deploys the site only for `v*` tags, so choose a name that cannot be mistaken for a release: **check that the Pages workflow would not fire on it**).
4. **Determinism check (phase 2):** yes or skip? Recommendation: yes, 200 genomes, about two hours.
5. **Homothallic policy outside Dothideomycetes:** apply as written, or restrict the first run to families where it was measured (Ascomycota MAT)? Recommendation: apply as written but report the count per phylum in the pilot.
6. **If the *Cochliobolus* scoped roster does not validate** in time: run without it and re-run only the Pleosporaceae, or wait? Recommendation: run without it; a Pleosporaceae re-run is small.

## 10. Risks
Queue contention (shared account); timeouts on large genomes (the 14400 s class); empty or missing FASTA (for example `GCA_986281625.1_DAH1041FM` has an empty file) listed rather than retried; a database or code change during the run (the frozen worktree prevents it); the intermittent test (`tests/detect/test_pipeline.py::test_fast_path_missing_gene_is_polished_in_the_existing_clusters_window`) also fails on unmodified main and is unrelated to the run; and the unverified items above.

## 11. After the run
Re-derive the family counts in `docs/notes/2026-10-18_families-needing-curation.md` from the new results (they come from v0.6.0); run the 12-locus sequence check for the strong-single-gene rule (decided to follow the clean run); re-run the Basidiomycota scope-only family table; update the failure tracker's revisit list with the new numbers.
