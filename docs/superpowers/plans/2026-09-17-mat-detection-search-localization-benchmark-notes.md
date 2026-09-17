# Task 12 benchmark/timing evidence: localize-then-polish search revision

Date: 2026-09-17. Environment: UCR HPCC interactive SLURM job (host `r11`),
pixi `default` environment with real `tblastn` (2.17.0), `makeblastdb`
(2.17.0), `miniprot`, `exonerate`, and `diamond` binaries on `PATH`
(confirmed via `pixi run which tblastn miniprot exonerate makeblastdb`).

This is a manual, one-time evidence log, not a pytest suite. Numbers below
are real measurements from real runs against a real genome assembly. Where
something was not measurable in this environment, that is stated plainly
rather than estimated.

## Step 1: chosen `assembly`-type curated record

`db/Basidiomycota/Agaricales/5334_h4-8_Aalpha_4/metadata.yaml` has
`locus.core.segments[0].sequence_source.type: assembly`:

- `accession: GCF_000143185.2` (*Schizophyllum commune* H4-8, RefSeq)
- `seq_region: NW_026089539.1`
- taxid `5334`
- curated locus: `NW_026089539.1:1821008-1827457` (genes `Z` 1821008-1824097(-),
  `Y` 1824369-1827457(+))

This assembly is small (25 scaffolds, ~39 Mb uncompressed FASTA, 12 MB
gzipped from NCBI's RefSeq FTP), which made a full real-genome run feasible
within a few minutes -- I did not need to fall back to a partial/synthetic
contig.

The genomic FASTA was downloaded directly from
`https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/143/185/GCF_000143185.2_Schco3/GCF_000143185.2_Schco3_genomic.fna.gz`
into node-local scratch (not committed to the repo -- it is a large,
re-downloadable third-party file, consistent with existing repo practice of
not vendoring full assemblies).

## Step 2: real before/after `matpredict detect` run

**Post-revision (current `main`, commit `8049a26`)**, real binaries, full
genome, `--taxid 5334`:

```
time pixi run matpredict detect --genome GCF_000143185.2_Schco3_genomic.fna --taxid 5334 --out-dir out_post
```

Result: completed in **1m36.4s wall clock** (1m23.8s user, 10.0s sys).
Output: `detected 13 candidate locus/loci (4 families attempted)`, with
`Basidiomycota:Abeta`, `Basidiomycota:Balpha`, `Basidiomycota:Bbeta`
honestly reported as not-detected ("no reference-protein hits found for
this family in this genome" -- correct, since this record only curates the
Aalpha sublocus for this species; Abeta/Balpha/Bbeta records exist for other
species, not this one).

Among the 13 Aalpha candidate calls, one matches the curated locus almost
exactly: `NW_026089539.1:1821007-1827453` (curated: `1821008-1827457`) with
gene `Z` polished by `miniprot_refine` at 98.31% identity and gene `Y` at
`tblastn_genome` (unpolished) 98.771% identity, both essentially the curated
gene models. Confidence tier reported: **medium**, not high -- capped
because gene `Y` never received a confirmed polish model from either
`exonerate --refine` or `miniprot` in this run (`status: unpolished`), which
is a real, reportable outcome of the tier logic (`_any_gene_unpolished`
caps at Medium regardless of raw identity), not a bug in this evidence run.
The remaining 12 candidate calls are lower-identity/lower-confidence
paralog-like hits elsewhere in the genome (expected: MAT pheromone-receptor
and HD genes have gene-family relatives outside the locus).

**Pre-revision (commit `cd7e386`, the last commit before Task 1 started the
localize-then-polish work)**, same genome, same taxid, same real binaries
(reused from the current pixi env -- the pre-revision code only needs
`exonerate`/`diamond`, which were already present; no separate env rebuild
was needed):

```
cd <worktree at cd7e386>
PYTHONPATH=<worktree>/src <pixi python> -m MATPredict detect --genome GCF_000143185.2_Schco3_genomic.fna --taxid 5334 --out-dir out_pre
```

(Run via a temporary `git worktree` at `cd7e386` plus `PYTHONPATH`, reusing
the main checkout's already-installed pixi binaries, rather than solving a
second pixi environment -- this kept the pre-revision run cheap to set up.)

Result: **did not complete within a 5-minute (300s) wall-clock budget** --
killed by `timeout 300` (exit 124/143), and no `detection_report.yaml` was
written. `ps` during the run showed the pre-revision pipeline had reached
its genomic path (`search_genomic`, `exonerate --model protein2genome`
against the *whole* 39 Mb genome with the *whole* curated reference-protein
set as query) and was still running exonerate alone at >3 minutes of CPU
time when the timeout fired. This matches the pre-revision code path read
directly from `cd7e386`'s `pipeline.py`/`search.py`: with no `--proteins`
supplied, the old code's only genome-only path was a single whole-genome
`exonerate --model protein2genome` splice-aware alignment of the entire
reference-protein FASTA against the entire genome FASTA -- there was no
tblastn localization step to cheaply narrow the search first.

**Timing comparison:**

| | wall clock | outcome |
|---|---|---|
| pre-revision (`cd7e386`) | **>300s (timed out, no result)** | no report produced |
| post-revision (`8049a26`, current `main`) | **1m36.4s** | 13 candidate calls, 1 near-exact match to curated Aalpha locus at medium confidence |

This is a real, if informal, confirmation of the localize-then-polish
revision's stated purpose: replacing an unbounded whole-genome
`exonerate --model protein2genome` scan (pre-revision) with a single batched
`tblastn` localization pass followed by narrow-window polishing
(post-revision) turned a genome-only run that did not finish in 5 minutes on
this ~39 Mb genome into one that finished in under 2 minutes. I did not let
the pre-revision run continue past 5 minutes to find its true completion
time -- per the task's own timeout guidance, forcing it to completion
was out of scope, and "it did not finish in the same budget the new code
finishes well within" is itself the load-bearing finding here, not a precise
pre-revision wall-clock number.

## Step 2b: what was NOT measured, and why

- **Exact pre-revision wall-clock time.** Only a lower bound is known
  (>300s). Getting an exact number would mean either running exonerate to
  completion (potentially tens of minutes on a 39 Mb genome against the
  order.yml-wide reference set for this taxid -- the whole reason for this
  revision) or profiling exonerate's own progress output, which it does not
  provide in a way `subprocess.run` output can be checked mid-flight without
  polling. Given the task's explicit few-minutes budget, I did not pursue
  this further.
- **A true recall/boundary-accuracy comparison via `matpredict detect
  benchmark`.** I read `src/MATPredict/detect/benchmark.py` in full. As
  currently implemented, `run_benchmark` only computes species/genus
  leave-one-out **holdout groupings** and reports `sensitivity=None` for
  every family in the curated DB -- the module's own docstring and inline
  comments say real per-species recall scoring against `run_pipeline` is
  "pending pipeline reference-injection support" and is an explicit,
  documented follow-up, not yet wired. I confirmed this by running it for
  real:

  ```
  pixi run matpredict detect benchmark
  ```

  Every family line printed `sensitivity=n/a`, either because there are
  `<=2` curated species for it, or because of the not-yet-wired recall
  scoring for families with more. **There is currently no real gene-recall
  or boundary-accuracy number the benchmark subcommand can produce for any
  family, pre- or post-revision** -- this is a genuine gap in the tool as it
  exists on `main` today, not something this task could work around, and not
  something to fabricate a number for. The spec's acceptance criterion (a)
  ("run the existing leave-one-out benchmark... to compare gene
  recall/boundary accuracy against a pre-revision baseline") could not be
  satisfied because the benchmark that would produce that comparison does
  not yet compute recall at all, in either the pre- or post-revision code
  (the benchmark module itself was not touched by this plan's 11 code
  tasks).
- **A second/third assembly-type record for cross-validation.** Only one
  `assembly`-type curated record exists in `db/` today
  (`5334_h4-8_Aalpha_4`); the other `type: assembly` metadata.yaml hits
  found by grep (`5346_a43-b43-okayama-7_HD_A43`, `5334_h4-8_Bbeta_2`,
  `5334_h4-8_Balpha_3`, `5346_a43-b43-okayama-7_PR_B43`,
  `5334_h4-8_Aalpha_4`) are all against the *same two* underlying genome
  assemblies (Schizophyllum commune H4-8 and Coprinopsis cinerea
  A43mut B43mut okayama7#130), so a genuinely independent second-assembly
  timing run was not attempted given the time budget for this task; the one
  real run above is the evidence this task produced.

## Step 3: honest summary

- Real, measured: post-revision `matpredict detect` against a real, full,
  ~39 Mb published genome assembly (GCF_000143185.2) with real `tblastn`/
  `miniprot`/`exonerate` binaries completes in 1m36s and recovers the
  curated Aalpha locus's gene coordinates almost exactly (medium confidence,
  due to one gene never reaching a confirmed polish model in this run --
  itself a real, reportable outcome, not a benchmark artifact).
- Real, measured: the pre-revision code's equivalent genome-only run against
  the same real genome did not finish in 5 minutes, consistent with reading
  its source: it had no genome-wide localization step and instead ran an
  unbounded whole-genome `exonerate --model protein2genome` scan.
- Not measured, and explicitly not fabricated: an exact pre-revision
  wall-clock number, and any gene-recall/Sn-Sp number from
  `matpredict detect benchmark` (pre- or post-revision) -- the benchmark
  subcommand does not yet compute recall for any family; that is a real,
  pre-existing gap in the tool, documented in `benchmark.py` itself as
  future work, not something introduced or fixable within this task.
