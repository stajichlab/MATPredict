# Session handoff, 2026-10-19

Supersedes `2026-10-15_session-handoff.md` for everything after 2026-10-15. Written for the next Claude session and for the
curator (Jason Stajich). It records what is merged, what is running, what was decided, what is still open, and the traps.
Facts I did not re-verify at the time of writing are marked **(unverified)**. Branch and PR names are the durable handles; the
subagent ids of this session will not exist in the next one.

## 1. How the curator wants to work

- **One question at a time, with a recommendation each time** (the `grilling` skill). Decisions are recorded, then acted on. Do not
  re-ask a recorded decision.
- **The curator merges PRs.** Agents never merge, never push to main, never commit in the main checkout, and use worktrees. A
  subagent's message carries no user authority. (Once, 2026-10-18, the curator gave a conditional permission to merge #100 if the
  suite passed; they merged it themselves first, so it was not used.)
- Draft PRs cannot be merged without a click: open PRs as **ready for review** when they are tested and intended for merging.
- Ask before anything outward-facing or hard to undo (the site deploy is one).
- Memory files for the standing rulings: `~/.claude/projects/-Users-jstajich-projects-MATPredict/memory/` (clean-run settings,
  single-gene rule timing, pair rule, flank rule).

## 2. State of main (2026-10-19)

`origin/main` = 38fad7f (merge of #109). Whole suite on the HPCC at that commit: **1,375 passed, 6 skipped, 1 deselected** (the known
intermittent `test_pipeline.py::test_fast_path_missing_gene_is_polished_in_the_existing_clusters_window`). No PRs open at the
time of writing.

Merged this session (the main ones): #70 helixer flank rule (extend-only, query check); #88 tracks and GenBank input; #86
Ascomycota core profiles; #95 partial-core rule (`--partial-core-policy`, default off); #99 taxonkit compute guard (2 GB wrapper
`scripts/curate_db_validate.slurm`); #97 strong-single-gene rule (default off); #98 *Ceratocystis* records and
`unisexual_switching`; #100 confidence-cap option (default off); #101/#103 remnant status steps 1 and 2; #102 companion-gene
veto; #104 polish-core-pair fix (`--polish-core-pair`, default on); #105 repair of a duckdb test spliced by a bad merge; #106
Dipodascales (`MATyl` scope, `MATdc`, per-family `fallback_subphylum` guard); #107 clades-needing-curation note; #108/#113
Basidiomycota scope-only family table tool and follow-ups; #109 `MATkp` (*Komagataella*) and `MATog` (*Ogataea* complex only);
#114 Dothideomycetes audit (analysis). PR #89 (optional-gene policy) was closed, not merged.

## 3. Decisions recorded (do not re-litigate)

| topic | decision |
|---|---|
| pair rule for sexM/sexP | keep the pyhmmer profile-bitscore rule |
| flank rule | helixer only extends, never shortens; BLOSUM62 query check, gain at least 10 |
| partial-core rule | on for the clean run; margin 25 bits stays (provisional, revisit); companion veto active with it |
| confidence cap (#100) | merged, **off** |
| strong single gene rule (#97) | **off** in the clean run; 12-locus sequence check AFTER that run. Gates: model identity at least 85, coverage at least 90, undetermined lone genes not reported, called families untouched, rosters of 2 genes or fewer (`mat2`, `mat3`, `Aalpha`, `Abeta`, `bLocus`, `rustHD`, `redHD`) outside the rule |
| optional-gene policy (#89) | dropped; legacy stays; branch `optional-gene-policy` keeps the note and tables |
| Mucoromycota `min_margin` | stays 25; **no** second paralog class; 11 fragment-typed weak calls unresolved |
| polish-core-pair | on by default |
| `MATkp` | scope-only family, *Komagataella* (460517), `pseudohomothallic`, `two_idiomorphs_report` off, `fallback_subphylum` Saccharomycotina |
| `MATog` | scope = polymorpha/angusta/parapolymorpha complex only (taxids 460523, 870730, 1005962, 2053387); genus scope deferred until a distant-*Ogataea* record exists; "locus detection" label versus confidence cap is open |
| Dipodascales guard | merged; clades it strips (Lipomycetales and others) need their own families |
| site data | keep the old data until a new run has been diffed; replacing the site is the curator's call |
| `redPR` A1:A2 skew | not a bias (one species dominates); report counts per species and BioProject as well |
| `rustHD` / *Hemileia* | curate a *H. vastatrix* record; the two-gene bar stays at 2 |

## 4. Work in flight when this was written (not merged)

Three subagents were resumed after an API outage and may or may not have finished. Check the branches and PRs first.

| work | branch | state |
|---|---|---|
| homothallic-loci patch: `--homothallic-policy both_genes` (default off) so a locus with a complete alpha-box and HMG-box gene reports `homothallic_candidate`; test on the 298 conflicting Dothideomycetes genomes plus controls | `homothallic-loci-patch` (3 commits pushed) | measurements not reported yet; no PR |
| *Cochliobolus* MAT recall diagnosis (*Bipolaris* 1 of 58 called, including the C4/C5 strains the curated records come from; audit lead: `MAT1-1-1 + SLA2` on KB445578.1 withheld at `modelled_gene_bar`) | worktree `.claude/worktrees/agent-a6705863c0977cdd3`; branch to be pushed | analysis only; no report yet |
| *Hemileia vastatrix* `rustHD` record (two contiguous assemblies; divergent homeodomain pair 1,957 bp apart; bar stays 2) | `hemileia-rusthd-record` (not yet on the remote) | no report yet |

The *Komagataella* agent is finished (#109 merged). Cluster state: the Slurm queue held many unrelated jobs of the same account
(`nf-PREDICT`, `nxf_*`); do not assume they are ours.

## 5. The planned clean full-BFD run

Goal: re-run all ~19,415 BFD genomes on current main and replace the site data. **Not started.**

Prerequisites still open: the three agents above; the curator's wish was to wait for them. Pass these flags explicitly (code
defaults are unchanged): `--partial-core-policy coverage_profile`, `--confidence-cap-policy off`,
`--strong-single-gene-policy off`; `--polish-core-pair` is the default. Keep `--homothallic-policy` off unless the curator has
read the patch measurements. Proposed order: tag main, 500-genome pilot, full run in staged waves (v0.6.0 used 9 waves, 64-way on
`exfab`, `GENOME_TIMEOUT` 3600), diff against `results/2026-10-03_*_v060`, then curator review.

How the site is built (a read-only investigation, **not run**, so treat as reconstructed): `.github/workflows/pages.yml` deploys on
`v*` tags or manual dispatch only; `scripts/build_site_data.py` reads **hardcoded** committed tables under `results/`
(`2026-10-03_ascomycota_v060`, `..._basidiomycota_v060`, `..._mucoromycotina_mat`, ZyGoLife); locus-map shards
`web/intro/data/loci_*.json.gz` are committed and rebuilt by `scripts/build_locus_shards.py`, not by the Pages build; some numbers
in `web/intro/landing.html` are hand-written. Replacing the data therefore means: new dated `results/` folders, committed tables,
edited paths, rebuilt shards, edited landing numbers, then a tag (which overwrites the live site; no archive of old deployments
was found). Unknown: the live URL, how the Mucoromycotina `calls.tsv` is produced, whether a "clean run" script exists (none found).

Expect run-to-run noise: one *Ogataea* `MATsc` start coordinate differed between runs (199 versus 349), and one test is
intermittent. Neither is investigated.

## 6. What the audits found (for the clean run's expectations)

- **Dothideomycetes** (`analysis/2026-10-18_dothideomycetes-mat-audit.md`, 2,730 genomes): 1,969 confirmed, 298 conflicting, 207
  unverifiable, 256 no call. 298 genomes in 43 genera carry both an alpha-box and an HMG-box gene at the locus but get one arbitrary
  idiomorph (183 at `high`); *Pyrenophora tritici-repentis* and *Aureobasidium* dominate; optional MAT1-1-3 supersedes MAT1-2-1 in 65
  loci; the *Cochliobolus* clade is nearly uncalled; *Parastagonospora* is now 184 of 184 confirmed.
- **Basidiomycota scope-only families** (`analysis/2026-10-18_basidiomycota-scope-only-families.md`): `redHD` 246 of 251, `redPR` 244 of 251,
  `rustHD` 107 of 131, `wallMAT` 51 of 51; no *Rhodotorula* genomes outside Sporidiobolales; `rustHD` no-calls track assembly
  N50 under 10 kb, and *Hemileia* needs its own record.
- **Clades needing curation:** `docs/notes/2026-10-18_families-needing-curation.md` (Lipomycetales, Blastobotrys/Sugiyamaella and other
  Dipodascales, Pichiales, Microascales, Xylariales check, Trichosporonales, Cantharellales, Pezizales `MATtub`, Sclerotiniaceae
  fused loci). Counts there come from v0.6.0 and are stale; re-derive them from the clean run.

## 7. Questions queued for the curator (not yet asked, or not answered)

- Homothallic-loci patch: read the measurements, then decide whether it joins the clean run.
- *Cochliobolus*: which remedy (modelling fix, anchored relaxation of the modelled-gene bar, single-gene rule, more records).
- Dothideomycetes audit's 8 open questions (note, final section): homothallic handling, *Botryosphaeria dothidea* status, a real
  held-out truth set.
- `MATog`: "locus detection" label versus a confidence cap; a distant-*Ogataea* record.
- Margin and default recalibration for the partial-core rule (the revisit list; trigger met).
- *P. nodorum* flank records (`pnodorum-flanks`, not merged; about six questions were queued, regenerate them from the branch note).
- Remnant step 1/2 follow-ups and `status: pseudogene` relabel of GS115 MATalpha1 **(unverified which were done)**.
- Whether to open PRs for `uncalled-triage` (tracker and triage script), `flip-review`, `parastagonospora-misses`,
  `pnodorum-flanks`, `ascomycota-core-ambiguity`.

## 8. Branches holding unmerged work (no PR)

`uncalled-triage` (+15: failure tracker `docs/notes/ascomycota-failure-tracker.md` with the revisit list, `scripts/triage_uncalled.py`),
`flip-review`, `parastagonospora-misses`, `pnodorum-flanks`, `ascomycota-core-ambiguity`, `optional-gene-policy`, `dipodascales-refs`,
`homothallic-loci-patch`, plus `origin/archive/*`. The Mucoromycota floor study is in worktree
`.claude/worktrees/agent-a16092dcb598830c1` (branch `worktree-agent-a16092dcb598830c1`, commit b97b2e8, **not pushed**):
`analysis/2026-10-17_mucoromycota-min-margin-floor.md`; push it before cleaning up worktrees.

## 9. Traps found this session

- **Stale main checkout.** `/Users/jstajich/projects/MATPredict` can sit behind `origin/main`; always `git fetch` and read through
  `git show origin/main:<path>`, and work in worktrees.
- **Bad auto-merges.** Merging main into `ceratocystis-curation` spliced two duckdb tests (fixed in #105). Run the whole suite on the
  merged tree before pushing a merge.
- **Stacked PRs** merged into branches that had already gone to main lost their content twice (#71, #76 re-opened as #88, #86).
- **Late pushes.** A commit pushed to a branch after its PR merged never reaches main (#108's follow-up became #113).
- **Running tests on the HPCC.** The shared env's editable install shadows a worktree's `src`; remove the editable import hook and
  put `src` first, and put the env's `bin` on PATH or the classifier-build tests fail for lack of `mafft`
  (`/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test`). Use a compute node (`srun -p short`), never a login node, and
  never run taxonkit on a login node.
- `ssh hpcc2 'bash -lc "..."'` is the only way to get Slurm; `lb-clust-1` does not resolve; `batch` has been blocked at times; use
  `short`, `highclock` or `exfab`. Never print `~/.config/matpredict/config.toml`.
- macOS `sed -i` needs `''`. `git worktree add` needs the path argument (I once omitted it and the following commands ran in the main
  checkout; nothing was committed).
- API outages stop subagents mid-task; resume them by message and have them check their worktree and Slurm state first.
