# Aligner comparison for the Mucoromycota sexM/sexP classifier (2026-09-29)

## Question
Should the classifier build swap MAFFT for FAMSA or MUSCLE 5? Curator additions:
which strategy `mafft --auto` picks; MAFFT E-INS-i; a trimming arm; per-arm HMM
relative entropy and match states. Goal: the most accurate HMMs, not speed.

## Data and code version
- Build code: polish-scope-cuts 64971b6 (deterministic build), exported read-only to
  code/; only the alignment step is swapped (common.py patches classifier_build._mafft).
- Training set: the shipped Mucoromycota set, sexP 76, sexM 9. P1 not involved.
- Positives: 85 training proteins + 23 Zygo locus proteins (108), 30 genera,
  genus-held-out; fragments generated exactly as results/2026-09-28_validation_f3_f4
  (full, HMG box +-5, three 50-90 aa windows; seed 11).
- Negatives: db/Mucoromycota/classifiers/MAT/paralog_negatives.faa (189); gate =
  nearest-rank 95th percentile of their best score (the shipped rule).
- Tools: MAFFT 7.526, MUSCLE 5.3, FAMSA 2.4.1, ClipKIT 1.3.0, pyhmmer 0.12.3. All arms
  single-threaded, sorted input, builder seed 42.

## Method
Per arm (run_arm.py): 3 repeat alignments and HMM builds per gene (byte identity);
full HMMs -> match states, alignment columns, mean match relative entropy, gate;
genus LOO over 108 positives x fragment types. Trimming arms on mafft_auto and FAMSA:
ClipKIT kpic-gappy; gap50 (drop columns >50% gaps); hmgbox (PF00505 envelope +-10
columns).

## Results
`mafft --auto` selects L-INS-i for both genes, so mafft_auto == mafft_linsi.

| arm | det | gate | sexP M/cols | sexM M/cols | full wrong / called>=25 / worst | HMG box called / worst | windows called / worst | true full >= gate | Zygo min margin |
|---|---|---|---|---|---|---|---|---|---|
| mafft_auto (=L-INS-i) | yes | 100.3 | 293/445 | 210/328 | 0 / 106 / 21.8 | 104 / 12.1 | 300 / 2.4 | 96 | 48.0 |
| mafft_einsi | yes | 101.0 | 293/444 | 199/323 | 0 / 106 / 23.1 | 104 / 13.6 | 299 / 3.8 | 96 | 51.0 |
| muscle5 | yes | 99.4 | 288/471 | 190/352 | 0 / 107 / 21.6 | 103 / 14.9 | 298 / 5.8 | 96 | 47.6 |
| famsa | yes | 100.7 | 303/448 | 201/306 | 0 / 106 / 23.1 | 103 / 16.6 | 301 / 6.7 | 96 | 45.6 |
| mafft_auto + gap50 | yes | 99.3 | 285/301 | 210/210 | 0 / 106 / 22.6 | 103 / 11.8 | 298 / 2.1 | 96 | 47.2 |
| famsa + gap50 | yes | 100.1 | 292/301 | 201/202 | 0 / 106 / 23.0 | 103 / 15.7 | 300 / 5.8 | 96 | 49.3 |
| mafft_auto + hmgbox | yes | 90.0 | 89/98 | 77/94 | 0 / 96 / 6.7 | 100 / 7.9 | 279 / 3.3 | 88 | 15.2 |
| famsa + hmgbox | yes | 90.6 | 95/130 | 81/94 | 0 / 101 / 11.2 | 102 / 11.1 | 281 / 1.5 | 88 | 25.0 |
| mafft_auto + clipkit | yes | 46.9 | 280/345 | 129/133 | **16** / 87 / 6.7 | 86 / 7.3 | 256 / 1.2 | 85 | **-33.2** |
| famsa + clipkit | yes | 44.9 | 283/329 | 126/133 | **16** / 87 / 2.3 | 86 / 2.8 | 257 / 0.2 | 85 | **-18.5** |

"called" = margin >= 25 bits; n = 108 full, 108 HMG box, 324 windows.

- Every arm is deterministic (3/3 identical alignments and HMMs).
- Mean match relative entropy is 0.589-0.591 for all untrimmed and gap50 arms;
  higher only for hmgbox (0.60-0.73), whose shorter models score worse.
- Untrimmed aligners are within noise of each other on every metric.
- ClipKIT kpic-gappy is harmful: 16/108 full proteins typed wrong, Zygo min margin
  negative, gate falls to ~45 bits. HMG-box restriction lowers accuracy and the
  gate. gap50 matches untrimmed.
- Runtime (arms ran concurrently on 2 cores, so relative only): FAMSA ~60 s,
  MAFFT L-/E-INS-i ~520-560 s, MUSCLE 5 ~1,670 s per full arm.

## What changed in detection
Nothing. No repository file was edited.

## Limits
- One classifier (Mucoromycota). sexM has 9 training sequences.
- Runtimes are under contention.
- Bioconda: famsa 2.4.1, muscle 5.3, mafft 7.525, clipkit 2.14.0 exist, but the
  anaconda API listing showed famsa/muscle builds for linux-aarch64/osx only in the
  last three subdirs returned; linux-64 availability was not confirmed.

## Curator decisions
- Open: keep MAFFT (accuracy tied; already in pixi) or switch to FAMSA (same accuracy,
  ~9x faster). Do not adopt ClipKIT or HMG-box trimming.

## Files
common.py, run_arm.py, res_<arm>__<trim>.json, loo_<arm>__<trim>.tsv, log_*.txt,
code/ (read-only export of 64971b6).
