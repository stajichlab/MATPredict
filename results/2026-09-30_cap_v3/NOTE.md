# V3 polish-cap rank (strong-fragment clusters first) -- first real regression check

Date: 2026-10-01. Status: awaiting curator sign-off (branch `cap-v3`, not in PR #9).

## Question
Does ranking strong-fragment clusters first inside the per-family polish cap
(V3, curator's ruling 2026-09-30) change any existing call, and does it recover
M. pusillus NRRL A-13674?

## Data and code version
- Baseline: frozen worktree `run-ad3ea9b` (PR #9 head ad3ea9b).
- Candidate: frozen worktree `run-cf5cb1e` (branch cap-v3, commit cf5cb1e = ad3ea9b + V3).
- Regression panel: `testset/regression_panel.tsv` (30 Ascomycota, 33
  Basidiomycota, 80 Mucoromycota, 20 record sources) + Zygo 23 on scaffolds and
  contigs, via `scripts/run_regression_panel.sh` (SLURM 29307033-29307043).
- LCG Mucoromycotina held-out set: 621 genomes, `--phylum Mucoromycota`, three
  chunks per side (`run_lcg.slurm`, SLURM 29307044-29307049).
- Every genome has a report on both sides; no run logged a failure.

## Method
V3 (`strong_fragment_scores`, `rank_for_polish`, `select_polish_clusters` in
`src/MATPredict/detect/pipeline.py`): a cluster whose unmodelled tblastn HSP
fragments score at or above the family classifier's gate threshold, and that no
paralog class claims, ranks first within the cap. The cap size never grows;
families without a classifier keep the plain rank. Diffs with
`scripts/regression_check.py diff`.

## Results
| panel | genomes | loci changed | changes touching a call |
|---|---|---|---|
| Ascomycota | 30 | 0 | 0 |
| Basidiomycota | 33 | 0 | 0 |
| Mucoromycota | 79 | 12 | 0 |
| record sources | 20 | 0 | 0 |
| Zygo 23 (both inputs) | 46 | 0 | 0 |
| LCG Mucoromycotina | 621 | 13 | 1 |

- The only call change: **M. pusillus NRRL A-13674** (scaffold_161) goes from
  `withheld:modelled_gene_bar+polish_capped` to **called Plus/medium**, typed
  from a gene model, margin 90.5 (fragment 97.5 before), best score 200.1.
- Zygo 23: 23/23 locus and idiomorph on scaffolds and on contigs, both sides.
- The other 24 changes (12 Mucoromycota panel, 12 LCG) are withheld-only: a
  different cluster now gets polished inside the cap (Synrac1, the two
  GCA_02363031x/2x Mucor genomes, Mycotypha Mycafr1, and LCG equivalents), and
  each stays withheld (modelled-gene bar or fraction floor). None became a call.
- Record self-check (candidate, all panel runs): 23 called, 1 missed
  (498019_b11221_MTL_alpha, the known GCA/GCF contig-naming false miss), 1
  withheld (5011_unknown-1_PM_combined, modelled-gene bar, already known), 6 no
  report, 87 no assembly accession.
- Runtime: base and candidate job walls were within a few percent per group
  (e.g. Mucoromycota 17:58 vs 16:34; LCG chunks 46-73 min vs 48-76 min).

## Launcher (first end-to-end run)
- It ran cleanly: 10 panel jobs, a dependent diff job, summaries written.
- Fixed: the self-check scanned only the `record_sources` group, so records whose
  source genome sits in another group showed `no_report`. It now scans every
  candidate group (`--reports $OUT/cand/*/runs`).
- Gap, not changed (panel is curator-owned): three record source genomes are not
  in the panel although they are in BFD -- GCF_000143185.2 (Schizophyllum H4-8:
  Aalpha/Balpha/Bbeta records), GCA_016772295.1 (Coprinopsis okayama-7: HD and PR
  records), GCA_964255735.1 (fen-198). These are the 6 `no_report` records.

## Limits
- V3 can only act in families with a classifier (Mucoromycota now).
- The benefit rests on M. pusillus here plus the 19/20 stress-test recoveries
  (results/2026-09-30_cap_protection/NOTE.md).

## Curator decisions
- Sign off V3 for PR #9 on this regression summary (1 call gained, 0 lost or
  changed; Zygo 23/23).
- Whether to add the three missing record source genomes to the regression panel.

## Files
- Summary: `regression/diff/summary.md`; per panel
  `regression/diff/<panel>/regression_summary.md` and
  `regression_withheld_changes.md`; TSVs `regression_diff.tsv`.
- LCG: `lcg_diff/regression_summary.md`, `lcg_diff/regression_withheld_changes.md`.
- Self-check: `regression/diff/record_selfcall_candidate_all.tsv`.
- Jobs: `regression/jobs.tsv`, `jobs_lcg.tsv`; launcher `run_lcg.slurm`.
