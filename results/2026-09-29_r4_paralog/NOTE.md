# R4: a non-MAT paralog class (P1) in the Mucoromycota idiomorph classifier

Status: decided (curator ruling J. Stajich, 2026-09-29); implemented and measured.

## Question
The weak "sexM" calls near 80 bits in Mucor genomes are a non-MAT HMG gene
("P1") present in both mating types (results/2026-09-29_sexM_like_paralog/).
Does a third classifier class for P1 withhold exactly those calls, without
collateral loss?

## Data and code version
- Code: branch r4-paralog, commit 41bd471 (parent 7a98c55 = PR #9 head, with the
  per-build gate threshold 9a458de and two_idiomorphs). Frozen worktree
  .claude/worktrees/run-41bd471.
- Classifier: db/Mucoromycota/classifiers/MAT. MAT HMMs NOT rebuilt (rebuilds are
  not bit-reproducible); P1 added with `build_idiomorph_hmms.py --paralogs-only`.
- P1 source: db/Mucoromycota/classifiers/MAT/paralogs/P1.faa (ONE sequence,
  GCA_000697295.1 Mucor indicus B7402, BFD; in no held-out set) + P1.yaml
  (provenance, reason, evidence, ruling, limit).
- Inputs: Mucoromycota 293 (results/2026-09-26_early_diverging/lists/Mucoromycota.tsv);
  Zygo 23 both inputs; LCG Mucoromycotina 621 (results/2026-09-28_lcg_holdout/lists2);
  Mucor_Jena 64 scaffolds (results/2026-09-28_mucor_jena_holdout/inputs.tsv).
- SLURM jobs 29207238-29207242 (jobs.txt).

## Method
- A classifier may list `paralog_classes` in its manifest. Each is built by the
  build script from a curated `paralogs/<name>.faa` with a required
  `paralogs/<name>.yaml` (source, reason, evidence). A one-sequence class is
  built from the sequence itself (pyhmmer Builder.build), as the replay did.
- A verdict records paralog scores; `paralog_class` is set when the best
  paralog score beats the best MAT score by >= min_margin (25 bits, the same
  typing margin).
- The gate step withholds such a call with reason `paralog_class` (before the
  MAT-gene gate; split-locus calls included), and lists it in `suppressed_loci`
  with its scores. Withhold, not `undetermined`: an undetermined call is still
  reported as a MAT locus, and P1 was the only call in 5 genomes.
- Build check: 0/85 MAT training proteins classed P1; the P1 source itself is.
  No leave-one-out is possible with one sequence.

## Results
compare_output.txt, changes.tsv (every withheld and changed locus).

| Set | Genomes | Called before -> after | Reported calls withheld (paralog_class) | Calls gained |
|---|---|---|---|---|
| Mucoromycota 293 (vs 3f755db = 9a458de code) | 293 | 253 -> 253 | 2 | 0 |
| LCG Mucoromycotina (vs 076afe4 runs2) | 621 | 536 -> 536 | 18 | 2 |
| Jena scaffolds (vs 076afe4) | 64 | 61 -> 61 | 6 | 3 |

- Zygo 23: 23/23 locus and idiomorph on scaffolds and on contigs (zygo23_41bd471/score.txt).
- LCG clean strain labels (n=17; disputed, misidentified and leaked strains
  excluded): agree 10, disagree 3, uncalled 4, before and after. No change.
- Literature positives (Syzygites, Zygorhynchus, M. genevensis, M. azygosporus): unaffected.
- No genome lost its only call (0 rows).
- Withheld reported calls match the replay: the 17 predicted LCG P1 calls plus
  one more, Mucor subtilissimus NRRL 6226 (Minus/medium; P1 149.8 vs sexM 85.7),
  which keeps a Plus/high call (sexP 297.8). All 6 predicted Jena calls.
- In Mucoromycota 293: the training genome GCA_000697295.1 (circular, expected)
  and GCA_060309335.1 (the BFD copy of Jena CBS210_80), which keeps its Minus/high call.
- Removing P1 lets the relaxed pass report the real locus in 5 genomes (all
  medium, partial_locus, model-typed, strong margins):
  - M. indicus NRRL 13468 and 13081: Plus (sexP 300.0, 301.3), sexP+rnhA.
  - Jena CBS221_71 Minus (sexM 168.2), CBS223_63 Minus (156.4), CBS763_74 Plus
    (sexP 309.6). CBS763_74 is one of the Jena strains that carried a strong
    sexP-like annotated protein elsewhere.
- The P1 class also labels 20+ previously unreported "undetermined" HMG loci
  (Syncephalastrum, Pilaira, Phycomyces NRRL 1555, Circinella umbellata,
  Benjaminiella; P1 ~80-106 vs MAT ~54-67) as paralog_class instead of their
  earlier withholding reason. No reported call is involved.

## What changed in detection
Commit 41bd471 (branch r4-paralog), fast-forwarded to polish-scope-cuts:
src/MATPredict/detect/{classifier,classifier_build,mat_gene_gate,pipeline,report}.py,
scripts/build_idiomorph_hmms.py (--paralogs-only), db/Mucoromycota/classifiers/MAT/
{manifest.yaml,paralogs/P1.faa,P1.yaml,P1.hmm}, tests/detect/test_paralog_class.py
(13 tests). Suite: 930 passed.

## Limits
- P1 is trained on ONE sequence. It is not validated beyond these sets.
- The P1 HMM also scores other HMG loci (e.g. Syncephalastrum) at ~100 bits:
  it is not specific to P1 in the strict sense. Only unreported loci were
  affected here, but a genome whose true MAT gene scores closer to P1 than to
  sexM/sexP by >= 25 bits would lose its call.
- LCG and Jena baselines are 076afe4 (before the per-build gate threshold and
  two_idiomorphs); changes were attributed by withheld reason, and every lost
  call is a paralog_class withholding.
- The training genome's own P1 call is withheld circularly.

## Curator decisions
- Made: R4 (2026-09-29).
- Open: whether P1 needs more training copies (other Mucor hiemalis/indicus
  genomes, outside held-out sets) before other families adopt paralog classes.

## Files
NOTE.md, compare.py, compare_output.txt, changes.tsv, jobs.txt, run_lcg_chunk.slurm,
run_jena.slurm, Mucoromycota_41bd471/, lcg_runs/, jena/runs_scaffolds/,
zygo23_41bd471/.
