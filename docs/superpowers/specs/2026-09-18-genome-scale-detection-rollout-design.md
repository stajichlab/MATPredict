# Genome-Scale Detection Rollout — Design

## Context and goal

Sub-project 2's `matpredict detect` pipeline (localize-then-polish search)
is implemented, reviewed, and empirically validated on 2 hand-picked real
genomes. Sub-project 1's reference database now covers real, tiered MAT-locus
records across Eurotiales (*Aspergillus*), Onygenales (*Coccidioides*),
Hypocreales (*Fusarium*), Sordariales (*Neurospora*), Helotiales
(*Sclerotinia*, *Botrytis*), several lichen families, and Pezizales
(*Tuber*), among others already in place.

Neither piece has been run at scale. `detect` takes exactly one genome per
invocation; nothing fetches real target genomes or aggregates results
across many of them; `benchmark.py`'s recall/precision scoring is a stub
that always reports `sensitivity=None`.

The goal of this next phase, per the user's own earlier framing: **"run
these predictions for say Ascomycetes, on the Eurotiomyces (et all the
Aspergillus and Penicillium) and see what our ability to extract these loci
is."** Concretely: acquire real genome assemblies for the taxa now covered
by the curated DB, run the existing detection pipeline against all of them,
and produce an honest, aggregated picture of detection performance —
surfacing further curation gaps and pipeline bugs the same way the 2
hand-picked validation runs did earlier, but at a scale that actually tests
the tool rather than 2 cherry-picked cases.

This is explicitly an **evaluation/rollout** phase, not a rewrite of
sub-project 2. `run_pipeline`/`search_localize`/`polish_with_exonerate`/
`polish_with_miniprot`/`tiering.assign_tier` are all staying as they are
unless this rollout finds a real bug in them (likely — the 2-genome
validation already found 2 Critical bugs this way; expect more at scale,
and treat that as the point of doing this, not a failure of the pipeline).

## Scope for this plan

**In scope:**
1. Genome acquisition for a defined target taxon list (starting with
   Eurotiomycetes: *Aspergillus*, *Penicillium*, plus the already-curated
   *Coccidioides*/Onygenales genomes as an immediate sanity-check tier,
   since ground truth already exists for them).
2. Batch orchestration: run `detect` across every acquired genome,
   correctly sized for this project's HPCC SLURM environment (see Global
   Constraints).
3. Result aggregation: one consolidated report across all runs — tier
   distribution, per-family detection rate, contigs/genomes where nothing
   was found, flagged anomalies (a family expected present but reported
   absent; a family reported at unexpectedly low tier).
4. A ground-truth sanity pass: for every target genome that IS (or is a
   close conspecific/strain match to) an existing curated record's own
   source genome, compare the detection result against that record's own
   curated coordinates — this is the cheapest, most reliable way to
   actually exercise `benchmark.py`'s currently-stubbed sensitivity
   scoring, using real data instead of a synthetic holdout.
5. A triage pass over whatever real gaps/bugs the rollout surfaces (new
   curation gaps, pipeline misses) — scoped as a fast-follow, not blocking
   the rollout's own completion.

**Out of scope (explicitly, do not let this creep in):**
- Ab initio prediction (sub-project 3) — the pipeline's `GeneEvidence`
  already carries a forward-compatible hook for this; this rollout does not
  need to build it.
- The ML/homology-corpus and ESM2 embedding ideas — separate future
  direction (see `matpredict_esm2_future_direction.md` memory), not this
  phase.
- Any change to the curated-record schema or curation workflow itself,
  unless a rollout finding makes one strictly necessary (e.g. a
  `taxonomic_scope` gap that silently mis-routes a whole target genome).
- A general-purpose genome-download tool for arbitrary taxa — build exactly
  enough to acquire this rollout's target list, not a reusable service.

## Target taxon list for the first rollout (proposed, confirm before
building)

Priority order, cheapest-ground-truth-first:

1. **Onygenales sanity tier**: the exact or closely related genome
   assemblies backing the newly-accepted *C. immitis*/*C. posadasii*
   records (`db/Ascomycota/Onygenales/`) — smallest, fastest, real ground
   truth already exists per-record.
2. **Eurotiales**: real *Aspergillus* genome assemblies beyond
   *A. nidulans*/*A. fumigatus* (their own source genomes are curated —
   include those too, as more sanity checks), plus *Penicillium* species
   (no *Penicillium* MAT record is curated yet — this is a genuinely blind
   test of the pipeline against a related-but-uncurated genus, exactly the
   kind of result that's useful).
3. **Sordariales/Hypocreales**: additional *Neurospora*/*Fusarium* species
   beyond the curated strains.
4. Broader Eurotiomycetes/Ascomycota, contingent on 1-3's results and
   however many findings need triaging first.

## Architecture

### Component 1 — Genome acquisition

A small, single-purpose script (not a general tool) that, given a target
taxon list (species name or taxid), fetches real genome assemblies from
NCBI. Use the `datasets` CLI (NCBI Datasets, if already available in this
HPCC environment — check `pixi.toml`/module system first rather than
assuming) or `efetch`/`esummary` against the assembly database if `datasets`
isn't already a project dependency. Do not build a second, competing NCBI
client — check whether `src/MATPredict/db/ncbi_client.py`'s existing
`NcbiClient` can be extended for assembly-level lookups (it currently
raises/skips on `GCA_`/`GCF_` accessions — noted in `validate.py`'s own
comments — this rollout may be the first real consumer that needs this
fixed, which is in scope since it directly blocks acquisition).

Downloaded genome FASTA files are large, per-genome; store per this
project's global storage rule (`.zst` for internal pipeline intermediates,
compressed by default). Do not commit downloaded genomes to git — they
belong in a git-ignored directory (`genomes/` or similar, added to
`.gitignore`) or, if run entirely inside SLURM jobs, in `$SCRATCH` with a
copy-back-what's-needed policy for anything kept afterward (see Global
Constraints).

### Component 2 — Batch orchestration

A wrapper that runs `matpredict detect --genome <fasta> --taxid <id>
--out-dir <dir>` (or calls `run_pipeline` directly, if that avoids
redundant CLI-parsing/setup overhead per genome — implementer's judgment,
grounded in whichever the Task 12 benchmark-notes timing data from
sub-project 2 suggests is the real bottleneck) across every acquired
genome.

**HPCC job sizing is a real constraint, not a suggestion**: per this
project's own established guidance, size each SLURM job toward ~1-1.5 hours
of real runtime, not one job per genome. The sub-project 2 benchmark notes
recorded a real single-genome runtime (~1m36s for a ~39 Mb genome
post-revision) — use that as the seed estimate, and bulk enough genomes per
job to hit the target window, adjusting per real timing once a handful of
jobs have actually run, not by guessing further. Do not default to
one-job-per-genome without first estimating whether that under-shoots the
target window.

### Component 3 — Result aggregation

Every `detect` run already emits a GFF3 + YAML report
(`report.write_detection_gff3`/`write_detection_report`). This component
reads every run's YAML report across the whole batch and produces ONE
consolidated summary: total genomes attempted, families detected per
genome, tier distribution per family across all genomes, genomes where a
curated family from the SAME order/class was expected but nothing was
detected (a real anomaly worth investigating, not necessarily a bug — could
be a genuinely divergent or absent locus, same as the real Xylariales
finding from this session's literature curation).

### Component 4 — Ground-truth sanity scoring

For any target genome that matches (exactly, or closely enough that a
domain judgment call is needed — flag ambiguous cases rather than silently
picking one) an existing curated record's own source genome/accession,
compare the detection result's reported coordinates against that curated
record's own coordinates. This is the first real, non-stubbed use of
`benchmark.py`'s intended sensitivity scoring — wire it up here rather than
leaving it a permanent stub, using real self-consistency data (a genome
scored against its own curated ground truth) as the first real test case
before attempting the harder cross-strain/holdout scoring `benchmark.py`
was originally scaffolded for.

**Binding design principle carried over from the exon/intron model work**:
this comparison must score on **protein-sequence match**, not exact
exon/intron structural agreement — reuse the same "protein match is ground
truth" principle here rather than reintroducing a naive coordinate-only
comparison that the exon/intron plan's Task 5 audit specifically confirmed
does not exist anywhere yet. Do not build a structural comparator; build a
protein-sequence one.

## Global Constraints

- HPCC SLURM job sizing: ~1-1.5 hours real runtime per job (see Component 2
  and this project's `nextflow-hpcc`-adjacent global guidance).
- Compress large intermediate/output files by default (`.zst` for internal
  pipeline intermediates, `.gz` only where a downstream tool needs it
  natively).
- Never hardcode `/scratch/$USER/...`; use `$SCRATCH` with `:?` so a missing
  var fails loudly, per this project's standing HPCC rules.
- Never commit downloaded genome FASTA files to git.
- Reuse `NcbiClient` for any new NCBI-facing code; do not build a second
  client.
- Ground-truth scoring compares protein sequences, never raw
  exon/intron/coordinate structure, per the binding design principle
  established in `docs/superpowers/plans/2026-09-18-mat-gene-exon-intron-
  model.md`.
- Follow this session's established review discipline for every non-trivial
  change: subagent-driven-development, independent live-data
  re-verification at every review step, never self-accept a curation
  proposal.

## What this spec deliberately does NOT settle (flag to the user before or
early in the implementation plan)

1. **Exact target genome count for the first rollout.** "A handful for
   sanity-checking, then Aspergillus/Penicillium broadly" is directional,
   not a number. The implementer should propose a concrete first-batch
   count (e.g. 10-20 genomes) once genome acquisition is working, rather
   than the plan hardcoding a number now.
2. **Whether `datasets` CLI is actually available in this HPCC environment.**
   Check first; the acquisition component's exact implementation depends on
   the answer.
3. **How ambiguous a "matches a curated record's source genome" call is
   allowed to get** before Component 4 should flag it rather than decide it
   automatically — needs a concrete example or two once real data is in
   hand, not a rule written in the abstract now.
