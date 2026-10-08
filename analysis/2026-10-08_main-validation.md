# main (e9060da) validation: full test suite, four real genomes re-run, report rendering
Status: open (tests pass; two rendering items and one input-handling item need a decision)

## Question
After PRs #45 to #48 (reads-type, A. fumigatus, per-genome HTML/PDF report), does `main` pass its tests, reproduce earlier campaign calls on real genomes, how long does it take, and do real campaign reports render with the pinned WeasyPrint 69?

## Data and code version
- Code: `main` e9060da (merge of PR #48), frozen worktree `.claude/worktrees/run-e9060da`, its own pixi `test` environment (`pixi install -e test --frozen`; WeasyPrint 69.0, pyhmmer 0.12.3, pytest 9.1.1).
- Inputs: four genomes from the 2026-10-03 campaigns, re-run with `detect` from the BFD input library and compared with the saved campaign reports (`results/2026-10-03_basidiomycota_v060/reports_all.tar.zst`, `results/2026-10-03_mucoromycotina_mat/reports_all.tar.zst`).
  Choice was mine, for variety, not random: Serpula lacrymans (HD + PR calls), Sporisorium reilianum (a/b loci), Microbotryum lychnidis-dioicae (no call), Absidia cylindrospora (Mucoromycotina, `--phylum Mucoromycota` as in the campaign).
- Jobs: tests 29613904 (8 shards) and 29614543 (re-run of one file); render part 1 29613905 (old-format YAML to HTML and PDF), part 2 29614533 (fresh `detect`, then report). All on partition exfab.
- The top-level checkout was fast-forwarded to e9060da. 64 untracked files that collided with files now on `main` were moved to `.claude/pull-backup-2026-10-08/` (63 identical to `main`; one, `results/2026-10-06_pr_array_support/compare_calls.py`, was an older version that uses the pre-rename `array_*` column names).

## Method
1. Test suite: the 120 `test_*.py` files split round-robin into 8 shards, one array task each (`results/2026-10-08_full_test_suite/run.slurm`), JUnit XML per shard.
2. Real-report rendering, part 1: `matpredict report genome --run <dir> --pdf` on the four saved campaign YAMLs (old format, no `run` block).
3. Part 2: `detect` on the four genomes (decompressed into `$SCRATCH` first), timed; then `report genome --pdf` on the fresh run directory. `compare_runs.py` compares loci (contig, start, end, idiomorph, confidence, locus class, genes found) and wall time.
4. One PDF was rasterised with Ghostscript (60 and 110 dpi) and read by eye (page 1 and 2 of the Serpula report).

## Results
Tests (`results/2026-10-08_full_test_suite/shards.tsv`):
- 1,126 tests in 8 shards. First pass: 1,119 passed, 7 failed, 0 errors, 0 skipped. All 7 failures are in `tests/detect/test_classifier_build_deterministic.py` with `RuntimeError: mafft not found on PATH`. Cause: my job called the environment's `python` by absolute path, so the environment's `bin/` (which holds mafft 7.526) was not on `PATH`. Re-run of that file with `PATH` set: 8 of 8 pass. So all 1,126 tests pass; the 7 were a job set-up error, not a code defect.
- Wall time per shard 336 to 429 s including start-up (pytest time 139 to 322 s), so the whole suite took about 7 minutes on 8 tasks. The earlier estimate of about 40 minutes was serial. Earlier counts were about 1,075 tests; the extra are the 36 report tests and the reads tests.

Accuracy of a fresh run against the campaign (`compare_runs.tsv`):
| Genome | Loci (campaign / fresh) | Identical loci, confidence, class, genes | Wall s (campaign / fresh) |
|---|---|---|---|
| Serpula lacrymans | 3 / 3 | yes | 144 / 87 |
| Sporisorium reilianum | 2 / 2 | yes | 111 / 70 |
| Microbotryum lychnidis-dioicae (no call) | 0 / 0 | yes | 330 / 198 |
| Absidia cylindrospora | 1 / 1 | yes | 149 / 54 |
- 4 of 4 genomes give identical loci. This is four genomes, not an accuracy estimate; the regression panel (`docs/regression-check.md`) is the gate for call changes and I did not run it.
- The fresh runs were faster in all four. The comparison is not controlled: the campaign ran 64 tasks at a time on other nodes, the fresh runs ran 4 at a time on one node alongside the test shards, and the code differs (the tiered polish and later changes). I did not isolate the cause.

Report rendering (`results/2026-10-08_report_render/out/`, WeasyPrint 69.0):
- All four old-format campaign YAMLs render to HTML and PDF without error. The fresh runs also write `report.html` by default.
- Fresh reports: HTML 28 to 78 KB; PDFs 11 pages (Serpula), 13 (Sporisorium), 21 (Microbotryum, no call) and 4 (Absidia); 42 to 88 KB.
- Timing: the first render of the old YAMLs took 522 to 546 s each (HTML and PDF in one call; four at once, on a node that was also running the test shards). The same step on the fresh runs, later, took 2 to 16 s. The 522 s is unexplained (cold start of the new environment on the shared filesystem is my guess; not tested). The HTML step alone was not timed separately.
- Defect seen in the Serpula PDF, page 2: the gene-evidence table is wider than its card; the table header and the "Reference record" column extend past the card's right border ([screenshot](figures/report_pdf_table_overflow.png)). I looked at pages 1 and 2 of this one PDF only. I did not view the other PDFs, the HTML in a browser, or a WeasyPrint 70 render, so I do not know if it also shows in HTML or in 70.
- The "Check: Idiomorph undetermined" note, MIP1 as a flank and the three-locus summary in the Serpula report match the campaign YAML; I did not read the other reports' text.

Input handling: `detect --genome` on a `.fa.gz` file fails with "Gzipped files are not suitable for indexing, please use BGZF" (Basidiomycota genomes) and `'utf-8' codec can't decode byte 0x8b` (the Mucoromycota genome), after about two minutes. The campaign scripts decompress to `$SCRATCH` first, which I did in part 2. The message does not say "decompress the file" (the earlier note on unreadable genomes found the same kind of unclear input error).

## What changed in detection
Nothing. This note adds results and one figure; no source changes.

## Limits
- Four genomes, chosen by me. Calls were identical on all four; this says nothing about the other 22,000 genomes.
- Timing is not controlled (see above); the render time of 522 s is unexplained.
- Visual check covers two pages of one PDF. WeasyPrint 70 (used by the other agent) was not run.
- Full test suite ran on one node type (gpu12) only.

## Curator decisions
Open: (1) table overflow in the PDF: fix in `src/MATPredict/report` (the other agent's code) or accept. (2) A clear error for gzipped genome input, or let `detect` read `.gz`. (3) Whether to re-run all BFD and the pangenome sets now; see the answer below. (4) Disk: `detect` now writes `report.html` per genome (28 to 78 KB here); for 22,685 genomes that is about 0.6 to 1.8 GB by extrapolation from four reports (not measured).

## Files
- `results/2026-10-08_full_test_suite/` (`run.slurm`, `rerun_classifier.slurm`, `shards.tsv`, logs).
- `results/2026-10-08_report_render/` (`genomes.tsv`, `run.slurm`, `run2.slurm`, `render.sh`, `compare_runs.py`, `compare_runs.tsv`, `out/*_fresh.html`, `out/*_fresh.pdf`, timing files).
- `analysis/figures/report_pdf_table_overflow.png`.
