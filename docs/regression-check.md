# Pre-sign-off regression check

Curator ruling, J. Stajich, 2026-09-30: before any new record is signed off,
and before any classifier rebuild ships, a replay of existing calls must list
every call whose gene model, label, confidence or presence changed. The curator
then judges each change. A code rule cannot tell a better gene model from a
worse one, so the check does not judge; it makes every change visible.

Why: in results/2026-09-29_umbelopsis_merge/, a new *S. racemosum* reference
changed the gene model detect built for *Circinella minor* (classifier 86.4 vs
111.8 bits) with no change to the classifier, and new *Umbelopsis* flank
references pushed *M. pusillus* out of the per-family polish cap. Neither was
visible from call counts alone.

## When to run it

- A new or changed curated record (any phylum).
- A classifier rebuild, or a new paralog class.
- Any change to a family roster, scope, cluster gap or detection rule.

## How to run it

1. Make two frozen worktrees on /bigdata: the baseline (current PR #9 head) and
   the candidate (baseline + the change). Never edit either while jobs run.
2. Submit the panel for both, plus the diff:

   ```bash
   BASE_WT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-<base> \
   CAND_WT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-<cand> \
   OUT=/bigdata/stajichlab/jstajich/projects/MATPredict/results/<date>_regression_<name> \
   bash scripts/run_regression_panel.sh
   ```

   This runs `testset/regression_panel.tsv` (163 genomes: 80 Mucoromycota incl.
   all Umbelopsidaceae, Syncephalastraceae and Lichtheimiaceae; 33
   Basidiomycota; 30 Ascomycota, one per order; 20 genomes that are the source
   of a curated record) and Zygo 23 on both inputs, one `short` job per group
   per side, then a dependent diff job.
3. Or diff existing run directories directly:

   ```bash
   python scripts/regression_check.py diff --out OUT --title "change vs baseline" \
       --pair Mucoromycota BASE/runs CAND/runs --pair Zygo23 BASE_ZYGO CAND_ZYGO
   ```

## What it reports

Per genome, loci are matched by family + contig + span overlap. Each locus on
either side gets one row in `regression_diff.tsv`. Change types:

| change | meaning |
|---|---|
| call_lost / call_gained | a reported call appeared or disappeared; the other side's state (absent, or `withheld:<reason>`, with `+polish_capped` when the cluster was capped) is shown |
| idiomorph_changed, confidence_changed, locus_class_changed, verification_changed | label changes on a call present on both sides |
| core_model_changed | a core gene's exon coordinates changed (length and identity shown) |
| gene_set_changed, span_changed | the genes found or the locus span changed |
| classifier_input_changed | model vs hsp_fragment typing changed |
| classifier_shift | classifier margin moved by >= 5 bits (`--score-delta`); smaller shifts stay in the TSV as `margin_delta` |
| withheld_reason_changed | the locus is withheld on both sides for a different reason |
| genome_missing_in_* | a genome has a report on one side only |

`regression_summary.md` lists every change that touches a call (a call on
either side). `regression_withheld_changes.md` lists every change to loci
withheld on both sides. Nothing is aggregated away.

## Sign-off rule

Attach `diff/summary.md` and the per-panel `regression_summary.md` to the
record's or rebuild's analysis report. The curator reviews every
`call_lost`, `call_gained`, `idiomorph_changed`, `confidence_changed` and
`core_model_changed` line before sign-off.

Validation: results/2026-09-30_regression_check_validation/ (a2fe1b4 vs
8eff88e) surfaces all four known changes -- the *Circinella minor* losses
(best score 111.8 -> 86.4, withheld by the MAT-gene gate), *M. pusillus*
(`withheld:modelled_gene_bar+polish_capped`), the *U. nana* gain and the
*M. griseocyanus* CBS 116.08 downgrade (Minus -> undetermined, margin 26.3 ->
17.8).
