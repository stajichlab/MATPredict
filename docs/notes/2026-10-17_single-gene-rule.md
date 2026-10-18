# Strong single gene rule (2026-10-17)

Branch `single-gene-rule`, draft PR #97. Code: `src/MATPredict/detect/strong_single_gene.py`, hooks in
`pipeline.py`, flag `--strong-single-gene-policy {off,on}` (default `off`, recorded in `run.parameters`).
Follows `2026-10-17_parastagonospora-misses.md` (branch `parastagonospora-misses`); the curator's decisions of 2026-10-17 are implemented as given.

## 1. What the rule does

With the policy `on`, in a family that has NO report in a genome after the strict, relaxed and split-locus passes, a cluster is reported when
all of these hold:

| condition | value | where |
|---|---|---|
| core gene role | `core_MAT` | `evaluate_strong_single_gene` |
| gene is modelled | status `polished_agree/disagree/single`, or an annotation hit (`diamond_proteome`), as in `split_locus` | same |
| identity gate | `GeneEvidence.identity` of the modelled gene >= **85** (curator) | same |
| near full length | `GeneEvidence.coverage` >= **90** (percent of the reference protein length; None fails) | same |
| lone | the cluster's evidence holds no other gene, except a hit that lost an idiomorph cross-match to it (`lost_to`) | same |
| idiomorph | determined (not `undetermined`) | same |
| MAT-gene gate | `apply_mat_gene_gate` in full: the classifier must have scored the modelled protein (no roster flank to fall back on), and the paralog-class check applies | pipeline |

The call is `idiomorph_gene_only`, confidence `low`, with a report block `strong_single_gene`
(`core_gene, core_identity, core_coverage, cluster_identity, reference_record, thresholds`).
`low` needed no change in the tiering code: `tiering.py` and the report writers carry the confidence as a string and nothing couples it to a class;
`report/html.py` ranks `low` and now adds a one-line note. A unit test writes the report row and checks `low` / `idiomorph_gene_only`.

Two or more qualifying clusters of one family: the best per idiomorph is kept, so two OPPOSITE idiomorph loci are both reported, each `low`.
The genome's `two_idiomorphs` statement is then written for that family (the roster flag `two_idiomorphs_report` is only set for Mucoromycota:MAT;
the pipeline now also enables the statement for any family holding a strong-single-gene call). No choice is made between the two and
`homothallism` is supported only by the same-locus arrangement, as before.

### Which identity (assumption A1)
The identity gate uses the MODELLED core gene's identity (`GeneEvidence.identity`, the polished model against its reference protein), the same
quantity `split_locus` gates on. The earlier agent used the cluster (tblastn hit) identity. They differ: the report carries both.
In the recovered loci the model identity is above the hit identity in most (99.39 -> 100.0, 85.61 -> 85.71, 88.05 -> 89.76, 86.4 -> 90.03) and below it in
GCA_035078975.1 (94.52 -> 88.62) and GCF_000523435.1 (94.29 -> 91.53). Cases with hit identity >= 85 but model identity < 85 are not reported and are not
counted anywhere in these tables (not measured; see questions).

### Why it cannot disturb other results (assumption A2, and the bug found by the control)
The first version admitted the lone cluster to polishing like any admitted cluster. The control sample then showed it was not neutral: in
GCA_003271955.1 (already called) and GCA_000568365.1 (uncalled) the extra polished cluster went through the ORDINARY path and produced a
`MATALPHA1/MATALPHA2` locus at **high** confidence, unflagged. That version is not what is proposed. The final version:

* the extra clusters ("ssg-only": core hit >= gate identity, not admitted by the floor, fraction found below the fraction floor and not ambiguous,
  so the strict loop would have dropped them without building) are polished on their localized genes only (no rescue search, so the cluster's hits do not change),
  do not compete for or consume the polish cap, and are skipped by the strict loop, the split-locus pass and the below-floor list;
* only `_strong_single_gene_calls` reads them.

Consequence: a pair that `off` would have BUILT (fraction found >= the floor, e.g. a one-gene roster) is not eligible. That excludes families with
a roster of two genes or fewer. I did not measure how many families that is (it is not Ascomycota:MAT, which has more than 4 genes).
A unit test compares every result of `on` that the rule did not make against `off`.

### Interaction with the other two policies
`--partial-core-policy coverage_profile` (on main): tested together (unit test; the rule still reports one low `MAT1-2` call).
`--idiomorph-policy {legacy,core_first}` (origin/optional-gene-policy, PR #89): not on main, so the unit tests skip there. I merged both branches
in a scratch worktree (conflicts only textual: adjacent parameters and the `_relaxed_results` signature), resolved them, and the new test file,
`test_partial_core_with_idiomorph_policy.py` and `test_idiomorph_policy.py` all passed (28 tests, including both `legacy` and `core_first` with the rule on).
That merge was not kept. The rule touches admission and a late pass only; `idiomorph_policy` acts inside `_build` (called by my pass) and in the overlap resolution.

## 2. Tests

`tests/detect/test_strong_single_gene.py` (17 + 2 conditional): off unchanged; on recovers the Parastagonospora pattern (lone MAT1-2-1 at 99.4 with a weaker cross-hit in its
cluster); 56, 70, 84.9% not reported; gate inclusive at 85 and decided on the model identity; short model (60% of reference) not reported;
unmodelled hit and flanking-only hit not reported; opposite-idiomorph pair reports both, both `low`/`idiomorph_gene_only`, `two_idiomorphs` written, no homothallism;
duplicate idiomorph reports the best only; an already-called family is left alone; no result that the rule did not make differs from `off`; confidence/class in the written report;
coexistence with `--partial-core-policy` (and `--idiomorph-policy` when present); CLI flag default `off`; unit tests of the evaluator.
HPCC (pixi env `test`, job sgr-tests, 2 CPUs, 16 GB): `tests/detect tests/db`: **1132 passed, 6 skipped** on the first version; final version: **1133 passed, 6 skipped** (the 6 include the 2 conditional tests above).

## 3. Evaluation

Method: `scripts/sgr_eval_driver.py` runs `matpredict detect` (lineage routing from the taxid, default parameters) twice per genome, policy `off`
then `on`, in one process. Search and polish calls are memoised between the two runs (deterministic functions of their arguments), so `on` is what an independent run would give
but costs a few seconds. Final `on` runs were re-run after the fix above with `off` kept from the first pass (`--on-only`).
Code: this branch on top of main 5b9a9e0 (v0.6.3). Cluster: 20 tasks, 2 CPUs, 8 GB each (the group memory limit held the 16 GB array pending), job names `sgr-*`.
Samplers: `scripts/sgr_make_lists.py` (seed 1). Raw rows: `/bigdata/stajichlab/jstajich/projects/sgr_run/out2/summary.*.jsonl` (not committed).

Note on "called/uncalled": the labels come from the 2026-10-03 v0.6.0 campaign; the code here is newer, so some "uncalled" genomes are called by `off`
(shown as "called by off").

### (a) The named genomes

| genome | `off` | `on` | model id | hit id | coverage |
|---|---|---|---|---|---|
| GCA_003497985.1 (P. nodorum) | none | MAT1-2 low `idiomorph_gene_only` | 100.0 | 99.39 | 100 |
| GCA_003503065.1 (P. nodorum, 318 bp from contig start) | none | MAT1-2 low | 100.0 | 99.39 | 100 |
| GCA_022559985.1 (P. avenae) | none | MAT1-1 low | 85.71 | 85.61 | 100 |
| GCA_003503165.1 (assembly gap) | none | none (as expected) | | | |
| GCA_035078975.1 | none | MAT1-1 low only | 88.62 | 94.52 | 99.0 |

GCA_035078975.1: the second idiomorph in the earlier note (MAT1-2 at 84.6%, hit identity) is below the gate, so the real genome gives ONE locus and no
`two_idiomorphs` statement. The opposite-pair behaviour is therefore demonstrated by the unit test only; no genome in this evaluation produced an opposite pair.
The gate at 85 is close to 84.6, so the outcome is sensitive to the identity definition (model vs hit) and to the threshold.

### (b) Pleosporales sample, 24 uncalled + 24 called (seed 1; NOT the earlier agent's list, which is not on the branch, so only GCA_035078975.1 and the three named genomes are certainly shared)

| | genomes | called by `off` | loci added by `on` | recovered | extra locus in called genome | unflagged changes |
|---|---|---|---|---|---|---|
| uncalled | 24 | 4 | 10 | 10 | 0 | 0 |
| called | 24 | 24 | 0 | 0 | 0 | 0 |

The 10 recovered loci: 8 MAT1-2, 2 MAT1-1; model identity 85.7-100 (median 91), coverage 97.4-100, all `low`.
The earlier offline filter reported 3 / 8 / 14 of 24 at hit-identity thresholds 95 / 90 / 80; 10 at 85 is in line with that.

### (c) Control, campaign-wide

| sample | genomes | called by `off` | loci added | genomes changed (off locus altered or lost) |
|---|---|---|---|---|
| called, 10 largest orders (60 each: Hypocreales, Eurotiales, Serinales, Saccharomycetales, Pleosporales, Botryosphaeriales, Pichiales, Glomerellales, Magnaporthales, Onygenales) | 600 | 600 | **0** | 0 |
| uncalled outside Pleosporales | 400 | 36 | 2 | 0 |
| BioSample truth set, `stated` MAT1-1 / MAT1-2 only | 203 | 200 | **0** | 0 |

Loci added in the uncalled sample: GCA_009802455.1 (MAT1-2, model 93.9, coverage 95.0) and GCA_011034965.1 (MAT1-2, model 91.9, coverage 100); neither has a truth label.
Truth set (203 genomes, MAT1-1/MAT1-2 labels only; the 58 MATa, 8 MATalpha and 27 homothallic/diploid rows were not run): `off` and `on` are identical: 178 agree, 21 disagree, 4 no call.
So the rule makes no call worse and none better on genomes with a truth label. **No added locus has a truth source**, so the correctness of the 12 recovered
loci (10 + 2, plus the three named) is unverified except by identity and, for the three Parastagonospora genomes, by the earlier agent's reading of the species.
Run time: `on` adds a few seconds per genome (cached search); mean `off` wall time per genome is about 80 s on these nodes.

Why the others are not recovered: not examined (identity under the gate, coverage under 90, an unmodelled gene, a `MAT-gene` gate withholding, or a family
the rule excludes). I did not tabulate it.

## 4. Assumptions and limits

1. A1 identity = modelled gene identity (above). A2 neutrality construction (above).
2. Coverage 90 is my choice (the curator asked me to state it); it equals the "full-length" bar of `partial_core`. The reported coverage is a LENGTH ratio, model length over the reference protein length,
   so a model with an internal insertion reads slightly high. The recovered loci have 95-100.
3. The idiomorph must be determined; an `undetermined` lone gene is not reported (my choice, to keep paralog-confusable calls out).
4. The full MAT-gene gate applies; a gated candidate goes to `suppressed_loci` with its gate reason.
5. Per idiomorph the best cluster is kept; a second copy of the same idiomorph is not reported.
6. Only families with no other report in the genome. A family already called is never extended, including with a second, opposite single-gene locus.
7. The split-locus limits and flank records are unchanged, as instructed.
8. The 60-per-order sample takes the whole order when smaller; random with seed 1 over the genomes of the campaign `genomes.tsv` that have an input FASTA and a taxid in the assembly summaries.
9. Sequence-level checks of the recovered genes (translation, domain check, synteny) were not done.
10. Cluster note: the memoised two-pass design gives identical results to two independent runs only if the polishing tools are deterministic; I did not test that directly.

## 5. Questions for the curator

See the PR description and the report; in short: model vs hit identity; 90 coverage; whether an `undetermined` lone gene should be reportable; whether a family that is already called should
accept a second, opposite single-gene locus; whether one-and-two-gene rosters need the rule.
