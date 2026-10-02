# Held-out rerun on current code, and curator species names (2026-10-01)

## Question
How do the two held-out Mucoromycotina sets score on the current PR #9 code, and
how should results be named? The curator supplied a species-name to CBS-number
table for the Jena strains.

## Data and code version
- Code: PR #9 `52b3ff9` (deterministic L-INS-i classifier, P1 paralog class,
  per-build MAT-gene gate, V3 cap ranking, two_idiomorphs, Umbelopsis and
  S. racemosum records). Frozen worktree `.claude/worktrees/run-52b3ff9`.
- Command as before: `matpredict detect --phylum Mucoromycota`, no taxid.
- Baseline: the 076afe4 runs (report 2026-09-28_heldout-sets.md).
- Held-out: no training or curation use.
- Curator table: `results/2026-09-28_mucor_jena_holdout/curator_names_cbs.tsv`.

## Method
- Curator names applied to both `curator_table.tsv` files (Jena file renamed from
  `per_strain.tsv`) by CBS, NRRL, RSA, IMI, URM, ATCC or BCRC number.
- Rerun via SLURM jobs 29318562-66; `analyze.py` for tables and label scoring;
  `scripts/regression_check.py diff` for every change.

## Results
### Curator names
- Jena: 59 of 64 strains named. Not in the table (taxonomy left blank):
  CBS169_57, CBS334_71, CBS417_77, CBS564_66, CBS608_78.
- LCG: 61 genomes renamed by collection number; 14 renames change the genus,
  e.g. Zygorhynchus -> Mucor (6), Mucor lamprosporus/recurvus/dispersus ->
  Backusella (5). LCG "Circinella rigida NRRL 2341" is Mucor durus, which
  explains why it was one of only two "Circinella" genomes with a full locus.
- Conventions: "T" = type strain; "T_of_X" = type strain of synonym X; NT, ET,
  IT, LT = neotype, epitype, isotype, lectotype.
- Ellisomyces anomalus NRRL 2465 Plus-T is the same strain as CBS 243.57 (type),
  putatively Plus. RSA_581- has an identical MAT locus (same coordinates, 307.8
  bits), likely the same strain or genome.
- Mating truth: none from the curator for either set, except Plus/Minus in file
  names (LCG 27, of which 7 are disputed). Jena folder names carry no labels.

### LCG (621 genomes)
| | 076afe4 | 52b3ff9 |
|---|---|---|
| Called | 536 | 536 |
| Plus / Minus | 247 / 246 | 257 / 249 |
| Both | 36 | 24 (all unlinked) |
| Undetermined | 7 | 6 |

- Best confidence (new): high 419, medium 115, low 2.
- Regression diff (`diff/lcg/regression_summary.md`): 22 call_lost (17 to the P1
  paralog class, mostly the M. indicus lineage incl. "B. ctenidia" NRRL 6239;
  5 to the MAT-gene gate, all Circinella incl. "R. microsporus" NRRL A-17693),
  5 call_gained (incl. M. indicus NRRL 13468 and 13081 Plus), 1 idiomorph
  change, 7 core_model_changed.
- File-name labels, clean set (no disputed, misidentified, training or Zygo):
  n=16, agree 11 (7 high, 4 medium), disagree 1, uncalled 4. The disagreement is
  Thamnostylum lucknowense RSA_1015_Plus-T (called Minus, medium). Uncalled:
  Fennellomyces x2, Thamnostylum x2. Before: 10 / 3 / 4.
- All 7 disputed labels still disagree. Zygo 23 subset: locus 23/23, idiomorph 23/23.

### Jena (64 strains)
| | 076afe4 | 52b3ff9 |
|---|---|---|
| Called | 61 | 61 |
| Plus / Minus | 26 / 31 | 29 / 31 |
| Both | 3 | 1 |
| Undetermined | 1 | 0 |

- Scaffolds and contigs give the same call: 63/63.
- Only CBS372_39 (Syzygites megalocarpus) carries both idiomorphs (unlinked).
- Changed calls: CBS763_74 (M. amphibiorum) Minus/high -> Plus/medium (real
  sexP locus, margin 254.5, after P1 withheld); CBS221_71 (M. ucrainicus) and
  CBS223_63 (Kirkomyces cordensis) Minus/high -> Minus/medium (real Minus loci,
  margins 102.4 and 89.4); CBS169_57 and CBS251_35 (M. falcatus) both -> Plus;
  CBS206_69 undetermined -> Minus/low.

### By species
`results/2026-10-01_heldout_rerun/lcg_by_species.tsv` (205 species; 3
misidentified genomes excluded) and `jena_by_species.tsv` (55 names).

## What changed in detection
None in this study (measurement only).

## Limits
- Median wall time ~260 s vs ~52 s, on unmatched nodes; a same-node runtime
  check is running (`results/2026-10-01_runtime_check/`).
- LCG names beyond the Jena-table matches cannot be verified; the curator has no
  framework to check them.
- Label truth is small (16 clean labels).

## Curator decisions
- Made (2026-10-01): rename to curator_table.tsv in both sets; apply curator
  names to LCG by collection number; keep file names elsewhere; mating truth
  only from file-name labels; rerun both sets.
- Open: the 5 unnamed Jena strains; T. lucknowense RSA_1015_Plus-T.

## Files
`results/2026-10-01_heldout_rerun/` (NOTE.md, summary.txt, lcg_per_genome.tsv,
lcg_by_species.tsv, jena_per_strain.tsv, jena_by_species.tsv, diff/);
`results/2026-09-28_mucor_jena_holdout/{curator_table.tsv,curator_names_cbs.tsv}`;
`results/2026-09-28_lcg_holdout/curator_table.tsv`.
