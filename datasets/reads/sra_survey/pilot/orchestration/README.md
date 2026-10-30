# Orchestration scripts of the raw-reads waves (archived 2026-10-28, issue #204)

These are the helper scripts that ran waves 1 and 2 and the v0.7.0 attribution work on the HPCC. They lived only in `/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit/` (the "scripts folder") and were not in the repository. The main pipeline pieces (`wave_driver.py`, `panel_from_assemblies.py`, `select_runs.py`, `genotype_mat_one_accession.sh`, `summarize_pilot.py`, `merge_rescues.py`, `clone_clusters.py`, ...) are one level up in `datasets/reads/sra_survey/pilot/`; the files here are the one-off drivers and probes around them.

**Read this before reusing anything.** They are archived as they ran: hard-coded HPCC paths (`/bigdata/stajichlab/...`, the frozen worktree `.claude/worktrees/run-b317ece`), no tests, many are short one-off probes. They were not re-run for this archive. The grouping below is by name and path, not verified script by script.

## Folders
- `wave2/`: wave 2 driving and finishing: `newpass.py` (species whose full stage was short of their cap), `mkfull3.py`, `mk_full2.py`, `fold3.py` (fold later passes into the final table), `wave2_final.py`, `wave2_final2.py` (final table with the all-strains *Botrytis cinerea* override), `wave2_clones.py` (clone pairs per species), `wave2_finish.sh`, `wave2_wait*.sh`, `yeastnames.py`, `w2tally.py`, `seven.py`, `mk_clonelist.py`, `submit_throttled.py`, `genus_lineage.py`, panel statistics (`mkpanel*.py`, `panelstats.py`, `pstats2.py`, `ssstats.py`, `po2.py`), `ogat.py`, `win.py`, and the *Botrytis* probes (`botr*.py`, `botcount.py`, `b05*.sh` for the B05.10 mapping test described in `../wave2/botrytis_B05.10_MAT_locus.md`).
- `wave1_and_pilot/`: small probes from wave 1 and the pilots (`asco5.py`, `ftp.py`, `rep.py`).
- `v070_run/`: the v0.7.0 attribution and integrity checks (`attr36*.sh/py`, `mk36.py`, `full_compare.py`, `seedcmp.py`, `r1r2b.py`, `chk1.py`, `integrity_v070.py`, `pack_v070.sh`, `tables_v070.slurm`, `dc.sh`, `dethelp.sh`, `run_asmcheck.sh`).
- `taphrinomycotina/`: the Taphrinomycotina candidate, polish and synteny work (`taph*.py`, `neo*.py`, `saito*.py`, `syn1.py`, `refs13.py`, `e2e*.`, `run_taph*.sh`, `dl.sh`); the maintained versions are `scripts/taph_mat_candidates.py`, `taph_polish_probe.py`, `taph_synteny.py`.
- `benchmarks/`: SRA throughput benchmarks and a compute-node probe (`bench*.sh`, `probe.sh`).

## Versions that differ from the repository
The HPCC scripts folder had older copies of three files than the repository: `clone_clusters.py` (52 lines differ; this is why the HPCC cluster job printed 219 clones for *B. cinerea* where the repository version gives 197 on the same pairs), `taph_mat_candidates.py` (30 lines) and `wave_driver.py` (2 lines: the `panel_no_unique_kmers` hold, added later). The repository versions are the maintained ones and are not duplicated here.

## Archive of the whole folder
The whole scripts folder (all scripts plus logs, small tables and run lists; reads, FASTA files, BAMs and sketches excluded) was packed on the HPCC before this commit:
`/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-28_v070_submit_scripts_archive.tar.gz` (10.6 MB), sha256 `732ec8830af90d53c244c8a1577fd849e20f4a238ecef846370ee981e150f4b1`. The folder itself was not touched or deleted. `MANIFEST_wave2.txt` records the run conditions. The wave 1 `waves/wave1_kexp/` result rows (k=21 and k=25 experiments) are only on the HPCC (`results/sra_reads_typing/waves/wave1_kexp/out`).

Nothing downloaded (reads, FASTQ, BAM) is committed.
