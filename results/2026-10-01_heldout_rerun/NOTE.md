# Held-out rerun on PR #9 52b3ff9 (Jena + LCG), by species under curator names

## Question
How do the two held-out Mucoromycotina sets score on the current PR #9 code
(deterministic L-INS-i classifier, P1 paralog class, per-build MAT-gene gate,
V3 cap ranking, two_idiomorphs, Umbelopsis/S. racemosum records), and what
changed since the 076afe4 runs?

## Data and code version
- Code: frozen worktree `.claude/worktrees/run-52b3ff9` (remote polish-scope-cuts
  52b3ff9). `matpredict detect --phylum Mucoromycota`, no taxid, as before.
- Jena: 64 genomes (scaffolds; 63 also on contigs). Curator names from
  `../2026-09-28_mucor_jena_holdout/curator_table.tsv` (no mating truth).
- LCG: 621 Mucoromycotina (`genomes/<Org>.sorted.fasta`). Curator names (61) and
  27 file-name labels (7 DISPUTED) from `../2026-09-28_lcg_holdout/curator_table.tsv`.
- Baseline: 076afe4 runs (`../2026-09-28_lcg_holdout/runs2`,
  `../2026-09-28_mucor_jena_holdout/runs_{scaffolds,contigs}`).
- SLURM jobs 29318562-66 (short); all completed; every genome has a report.
- Held-out: no training or curation use.

## Method
`run_lcg.slurm`, `run_jena.slurm` (detect); `analyze.py` (per-genome tables,
by-species summaries, label scoring, Zygo 23); `scripts/regression_check.py diff`
for the full change list (`diff/`).

## Results
LCG (621): called 536 -> 536. Call split old -> new: Plus 247 -> 257, Minus 246 ->
249, both 36 -> 24, undetermined 7 -> 6, uncalled 85 -> 85. Best confidence (new):
high 419, medium 115, low 2. two_idiomorphs: 24, all unlinked.
- Regression diff (call-touching): 22 call_lost, 5 call_gained, 1 idiomorph_changed,
  1 confidence_changed, 7 core_model_changed (`diff/lcg/regression_summary.md`).
  - Lost: 17 to `paralog_class` (P1) — M. indicus x8 (incl. B. ctenidia NRRL 6239,
    misidentified M. indicus lineage), M. hiemalis x2, M. racemosus, M. rouxianus,
    M. rouxii, M. luteus var. indica, M. subtilissimus, Mucor sp. x3; 5 to the
    MAT-gene gate — C. minor x2, C. umbellata NRRL 2417, R. microsporus NRRL A-17693
    (misidentified C. minor) — all Lichtheimiaceae/Circinella.
  - Gained: M. indicus NRRL 13468 and 13081 (Plus; real locus revealed once P1 is
    withheld), M. ramannianus NRRL A-21216 and Mucor sp. NRRL 1454 (Minus),
    S. racemosum NRRL 2495 (Minus).
  - M. pusillus NRRL A-13674: called Plus in both runs (V3 keeps it).
- File-label scoring, clean set (16; excludes 7 DISPUTED, misidentified,
  training-leak Pilaira RSA_1997_Plus, and Zygo): agree 11 (high 7, medium 4),
  disagree 1 (Thamnostylum lucknowense RSA_1015_Plus-T called Minus, medium),
  uncalled 4 (Fennellomyces x2, Thamnostylum x2). Labelled Zygo 3/3. DISPUTED
  7/7 still disagree (expected).
- Zygo 23 (LCG scaffolds): locus 23/23, idiomorph 23/23.
- By species: `lcg_by_species.tsv` (205 species; 3 misidentified genomes excluded).

Jena (64): called 61 -> 61. Split old -> new (scaffolds): Plus 26 -> 29, Minus
31 -> 31, both 3 -> 1, undetermined 1 -> 0, uncalled 3 -> 3. Scaffolds vs contigs
same call 63/63. two_idiomorphs: CBS372_39 (Syzygites megalocarpus), unlinked.
- Calls changed: CBS169_57 both -> Plus (2 P1 withheld); CBS251_35 (M. falcatus)
  both -> Plus (P1); CBS763_74 (M. amphibiorum) Minus/high -> Plus/medium (P1
  withheld; real locus sexP margin 254.5); CBS221_71 (M. ucrainicus) and CBS223_63
  (Kirkomyces cordensis) Minus/high -> Minus/medium (P1 withheld; real Minus
  locus, margins 102.4 and 89.4); CBS206_69 (Mucor sp. nov.) undetermined -> Minus/low.
- By species: `jena_by_species.tsv` (55 names; 5 strains not in the curator table).

Runtime: median ~260 s/genome on LCG (076afe4: ~52 s); L-INS-i classifier and
added steps; not node-matched, so indicative only.

## What changed in detection
Nothing changed here; this is a measurement of PR #9 52b3ff9.

## Limits
- LCG truth = file-name labels only (16 clean); Jena has no mating truth.
- Species names: curator names for 61 LCG and 59 Jena genomes; others are file names.
- Diff counts include withheld-only changes; call-touching changes are listed above.

## Curator decisions
Made: rerun both on current code (2026-10-01). Open: none raised here; the clean
disagreement (T. lucknowense RSA_1015_Plus-T) and 4 uncalled labelled genomes
remain for review.

## Files
`lcg_per_genome.tsv`, `lcg_by_species.tsv`, `jena_per_strain.tsv`,
`jena_by_species.tsv`, `summary.txt`, `diff/summary.md`,
`diff/{lcg,jena_scaffolds,jena_contigs}/regression_summary.md` and
`regression_diff.tsv`, `analyze.py`, `run_lcg.slurm`, `run_jena.slurm`, `jobs.txt`,
`lcg_runs/`, `jena_scaffolds/`, `jena_contigs/`.
