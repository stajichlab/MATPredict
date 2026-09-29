# Deterministic classifier builds
Status: decided (code); open (whether to adopt the candidate full rebuild)

## Question
Rebuilding the Mucoromycota classifier from the same training set moved
held-out scores by up to 8 bits (results/2026-09-29_gate_threshold/NOTE.md).
Where does that come from, can the build be made deterministic, and what would
a deterministic full rebuild change against the shipped HMMs?

## Data and code version
- Code: branch deterministic-build from polish-scope-cuts da35821; code commit
  a100cea.
- Database: the same tree. Shipped classifier db/Mucoromycota/classifiers/MAT/
  (built 2026-09-26 with the old code) is NOT changed.
- Replay inputs: results/2026-09-29_r4_paralog/Mucoromycota_41bd471/runs (293
  genomes, code 41bd471 = the shipped HMMs); genomes from BFD
  input_clean_genomes.

## Method
1. Tested MAFFT on the real 76-protein sexP training set: three runs with
   `--thread 4`, three with `--thread 1`, and one with input order reversed.
2. Fixed each source in src/MATPredict/detect/classifier_build.py:
   single-threaded MAFFT (`MAFFT_OPTIONS`), inputs sorted by id, a fixed
   pyhmmer Builder seed (42), a fixed HMM DATE (1970-01-01), no COM line; the
   manifest now records `deterministic`, `tool_versions` (pyhmmer 0.12.3,
   HMMER3/f [3.4 | Aug 2023], MAFFT v7.526), `mafft_options`,
   `hmm_builder_seed` and `training_sha256`.
3. Built the full Mucoromycota:MAT classifier twice, in two fresh processes,
   to candidate/ and candidate_b/.
4. Scored the shipped and candidate final HMMs on the 85 training proteins,
   the 112 Zygo 23 locus proteins (held out) and the 189 HMG paralog
   negatives (compare_candidate.py).
5. Replayed the 253 model-typed Mucoromycota calls: rebuilt each call's core
   proteins from the report's exon coordinates, scored with both classifiers,
   and compared verdicts (typing margin 25, P1 paralog class, MAT-gene gate
   with each build's threshold, >=2 roster flanks at >=40% as the flank route).

## Results
- MAFFT `--thread 4`: three runs, three different alignments. `--thread 1`:
  three runs, one alignment. Reversed input order: a different alignment.
  This is the cause of the rebuild drift.
- Two candidate builds in separate processes: byte-identical sexP.hmm,
  sexM.hmm and paralogs/P1.hmm (sha256 98f042b0..., a8d08271..., a802fe1a...),
  identical leave-one-genus-out scores and training checksum (ccfbf37a1ed5...).
  Each build took 3 min 23 s, 146 MB peak.
- Candidate vs shipped (compare_output.txt, score_deltas.tsv):

| Set | n | best-score delta (bits) | mean abs delta | sign flips | at/above gate (shipped / candidate) |
|---|---|---|---|---|---|
| training proteins | 85 | -8.7 to +7.6 | 3.30 | 0 | 80 / 80 |
| Zygo 23 locus proteins (held out) | 112 | -1.2 to +6.6 | 0.69 | 4 (non-MAT proteins at 0-1.2 bits) | 21 / 21 |
| HMG paralog negatives | 189 | -1.2 to +3.7 | 0.26 | 0 | 10 / 9 |

- Leave-one-genus-out: shipped 85/85, worst correct margin 22.5; candidate
  85/85, worst 21.8 (recommended min_margin 11.2 vs 10.9). Gate threshold
  from the build: shipped 99.9, candidate 100.3 bits.
- Replay: 253 model-typed calls; the shipped scores reproduced within 0.15
  bits for 177 (the other 76 differ because the pipeline scores both the
  miniprot and exonerate models while the report keeps one). Verdict changes
  candidate vs shipped: 0 of 177 reproduced, 0 of 76 on the proxy.

## What changed in detection
Commit a100cea: deterministic build code, build-script help on when to use a
full build, `--gate-only` and `--paralogs-only`, and 7 new tests
(tests/detect/test_classifier_build_deterministic.py). No shipped HMM changed.

## Limits
- The 76 unreproduced calls are compared on the one report model only.
- The replay covers Mucoromycota (293) only, not LCG or Jena.
- The Zygo file holds all 112 HMG-region proteins at the Zygo loci, not only
  the 23 MAT proteins.
- Determinism holds for the recorded tool versions; a different MAFFT or
  pyhmmer version may change the alignment or HMM.

## Curator decisions
- Made 2026-09-29: option (c) — deterministic builds plus the fast paths.
- Open: adopt the candidate full rebuild as the shipped classifier? It moves
  training-protein scores by up to 8.7 bits (mean 3.3) and the gate threshold
  from 99.9 to 100.3, but changes no verdict in the 253-call replay. Adopting
  it makes the shipped HMMs reproducible from the recorded inputs; not
  adopting it keeps the scores the signed-off records were measured on.

## Files
compare_candidate.py, compare_output.txt, score_deltas.tsv,
replay_changes.tsv (empty: no changes), candidate/ and candidate_b/ (candidate
HMMs and manifests; not shipped), candidate.log, candidate_b.log.
