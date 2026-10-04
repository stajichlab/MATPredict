# Full Basidiomycota run on v0.6.0 (2026-10-03)

All BFD Basidiomycota genomes with a taxid and a genome file (samples.csv), detect with
lineage routing from the taxid. Code: frozen worktree run-7c7ed99 (= v0.6.0).
Gate passed before launch: regression panel + Zygo 23, 0 changed loci
(../2026-10-03_regression_v060/diff/summary.md); pilot ../2026-10-03_basidiomycota_pilot.
Jobs (partition exfab): small_0..2 = genomes < 500 Mb, orders interleaved, 64-way,
GENOME_TIMEOUT 3600; large_0 = genomes >= 500 Mb, 32-way, GENOME_TIMEOUT 14400 (account exfab).
Each job writes <chunk>/runs, rollout_summary.yaml and reports.tar.zst.
Previous full run (older code run-ad1f865): ../2026-09-26_basidiomycota_full.
