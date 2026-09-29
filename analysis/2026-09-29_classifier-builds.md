# Classifier builds: aligner choice, determinism, gate threshold, P1 paralog class
Status: decided (aligner, determinism, gate, P1); running (deterministic full rebuild with explicit L-INS-i, awaiting curator approval)

## Question
How should the Mucoromycota sexM/sexP classifier HMMs be built so that they
are accurate and reproducible, and how should a non-MAT HMG paralog be kept out
of MAT calls?

## Data and code version
- Build code: `src/MATPredict/detect/classifier_build.py`,
  `scripts/build_idiomorph_hmms.py` on PR #9 (`polish-scope-cuts`).
- Shipped classifier: `db/Mucoromycota/classifiers/MAT/` (sexP n=76, sexM n=9
  training proteins; `paralog_negatives.faa`, 189 HMG paralogs;
  `paralogs/P1.faa` + `P1.yaml`).
- Positives for testing: 85 training proteins + 23 Zygo locus proteins (108),
  30 genera, genus held out.

## Method
1. Gate threshold per build (9a458de): each build scores the 189 paralog
   negatives and takes the nearest-rank 95th percentile as the MAT-gene gate;
   written to `manifest.yaml`; `--gate-only` computes it for existing HMMs.
2. Determinism (a100cea): find and remove sources of run-to-run change; prove
   with two builds in separate processes.
3. Aligner comparison (read-only): swap only the alignment step; compare MAFFT
   (`--auto`, L-INS-i, E-INS-i), MUSCLE 5, FAMSA, and trimming (gap50,
   HMG-box restriction, ClipKIT kpic-gappy); 3 repeats per arm.
4. P1 paralog class (41bd471): optional per-family paralog classes built from a
   curated source file with provenance; a call whose best class is the paralog
   by >= 25 bits is withheld (`paralog_class`); added with `--paralogs-only`.

## Results

### Aligner comparison (curator asked for this finding to be recorded)
Source: `results/2026-09-29_aligner_comparison/NOTE.md`.
- **No aligner is more accurate than another for these HMMs.** `mafft --auto`
  already selects L-INS-i (a high-accuracy mode) for both genes.
- **Trimming hurts.** ClipKIT kpic-gappy typed 16/108 full proteins wrong and
  gave negative Zygo margins; HMG-box restriction lowered accuracy.

| Arm | Gate (bits) | Full: wrong / called>=25 / worst margin | HMG box called / worst | Windows called / worst | Zygo min margin |
|---|---|---|---|---|---|
| MAFFT --auto (= L-INS-i) | 100.3 | 0 / 106 / 21.8 | 104 / 12.1 | 300 / 2.4 | 48.0 |
| MAFFT E-INS-i | 101.0 | 0 / 106 / 23.1 | 104 / 13.6 | 299 / 3.8 | 51.0 |
| MUSCLE 5 | 99.4 | 0 / 107 / 21.6 | 103 / 14.9 | 298 / 5.8 | 47.6 |
| FAMSA | 100.7 | 0 / 106 / 23.1 | 103 / 16.6 | 301 / 6.7 | 45.6 |
| MAFFT + gap50 | 99.3 | 0 / 106 / 22.6 | 103 / 11.8 | 298 / 2.1 | 47.2 |
| MAFFT + HMG-box trim | 90.0 | 0 / 96 / 6.7 | 100 / 7.9 | 279 / 3.3 | 15.2 |
| MAFFT + ClipKIT | 46.9 | **16** / 87 / 6.7 | 86 / 7.3 | 256 / 1.2 | **-33.2** |

n = 108 full proteins, 108 HMG boxes, 324 windows (50-90 aa). Every arm was
deterministic (3/3). Mean match relative entropy 0.589-0.591 for untrimmed and
gap50 arms. Runtime (relative only, shared cores): FAMSA ~60 s, MAFFT ~520-560 s,
MUSCLE 5 ~1,670 s.

### Determinism
Source: `results/2026-09-29_deterministic_build/NOTE.md`.
- Cause of the earlier up-to-8-bit rebuild drift: MAFFT `--thread 4` (three
  runs, three alignments), input order, and timestamps in the HMM files.
- Fix: single-threaded MAFFT on id-sorted input, builder seed 42, fixed HMM
  date, tool versions and a training-input checksum in the manifest.
- Two builds in separate processes: byte-identical sexP, sexM and P1 HMMs.
- Candidate full rebuild vs shipped: training scores -8.7 to +7.6 bits (mean
  |delta| 3.30); Zygo proteins -1.2 to +6.6 (0.69); paralogs -1.2 to +3.7
  (0.26); LOO 85/85 in both; gate 99.9 -> 100.3 bits; 0 verdict changes on 253
  Mucoromycota calls (177 reproduced exactly, 76 checked on a proxy).

### Gate threshold per build
Source: `results/2026-09-29_gate_threshold/NOTE.md`.

| Build | Threshold | Negatives at/above | Held-out positives at/above | LOO |
|---|---|---|---|---|
| PR #9 | 99.9 | 10/189 | 75/85 | 85/85 |
| curation-umbelopsis | 98.9 | 9/189 | 76/88 | 88/88 |

0 call changes on 293 Mucoromycota genomes in both trees; Zygo 23/23 on both
inputs.

### P1 paralog class (R4)
Source: `results/2026-09-29_r4_paralog/NOTE.md`.

| Set | Called before -> after | Withheld as paralog | Real loci revealed |
|---|---|---|---|
| Mucoromycota 293 | 253 -> 253 | 2 | 0 |
| LCG Mucoromycotina | 536 -> 536 | 18 | 2 |
| Jena (scaffolds) | 61 -> 61 | 6 | 3 |

- No genome lost its only call; Zygo 23/23 on both inputs; clean strain labels
  unchanged (10 agree, 3 disagree, 4 uncalled; n = 17).
- Revealed loci: M. indicus NRRL 13468 and 13081 Plus (sexP ~300); Jena
  CBS221_71 Minus (168.2), CBS223_63 Minus (156.4), CBS763_74 Plus (309.6).

## What changed in detection
- 9a458de gate threshold per build; a100cea deterministic builds; 41bd471 P1
  paralog class (all on PR #9). No shipped MAT HMM was rebuilt.
- Pending (branch `aligner-default`, `results/2026-09-29_aligner_default/`,
  running at the time of writing): explicit MAFFT L-INS-i
  (`--localpair --maxiterate 1000`) as the single default, then one
  deterministic full rebuild of sexM/sexP replayed on Mucoromycota, Zygo, LCG
  and Jena.

## Limits
- One classifier (Mucoromycota); sexM has 9 training sequences.
- P1 is trained on one sequence (GCA_000697295.1, M. indicus B7402); it also
  scores ~80-106 bits on 20+ other HMG loci, now labelled paralog_class.
- Aligner runtimes were measured under contention. Linux-64 bioconda builds
  for FAMSA and MUSCLE 5 were not confirmed.

## Curator decisions
- Made 2026-09-29: gate threshold from each build; deterministic builds with
  `--gate-only` / `--paralogs-only` fast paths (option c); build P1 class (R4);
  no full rebuild of a shipped classifier without approval; keep MAFFT with
  L-INS-i set explicitly; no trimming.
- Open: approve the deterministic full rebuild after its replay.

## Files
`results/2026-09-29_aligner_comparison/` (NOTE.md, res_*.json, loo_*.tsv);
`results/2026-09-29_deterministic_build/` (score_deltas.tsv, replay_changes.tsv);
`results/2026-09-29_gate_threshold/`; `results/2026-09-29_r4_paralog/`
(changes.tsv).
