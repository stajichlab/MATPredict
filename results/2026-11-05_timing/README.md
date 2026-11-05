# Timing benchmark, 20 genomes (2026-11-05, #258)

**Question:** the v0.8.0 definitive run took a median of 96 s per genome against 34 s for v0.7.0 on the same 19,343 Ascomycota genomes (ratio 3.9; 6.6 at the 90th percentile; `definitive_run_ratio.py`), and even genomes that were not searched at all were 7x slower. Is the v0.8 code slower?

**Method** (`run_benchmark.sh`, `compare.py`): 20 genomes (6 not_searched, 10 lineage-routed, 4 phylum-fallback; 9-51 Mb, 5-2,381 contigs) run with the v0.7.0 code and database (`--preset v0.7`) and with the v0.8.0 code and database (`--preset v0.8`) at the same time, 2 CPUs each, `wall_seconds` from the reports.

**Result** (`benchmark_20_genomes.txt`): median 19 s (v0.7.0) against 22 s (v0.8.0), median ratio 1.00; every genome within 0-2 s except two (+1-3 s on 2-9 s runs). **The code is not slower.** The 3.9x in the definitive run was cluster contention (shared storage and many jobs at once), not detection. A v0.8.1 full run should cost about what v0.7.0 did (151 h of per-genome wall time in total) on a quiet cluster, and up to about 560 h on a busy one.
