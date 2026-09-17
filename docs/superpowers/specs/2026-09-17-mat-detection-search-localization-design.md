# MATPredict Detection Pipeline: Search Stage Revision — Localize-then-Polish

Date: 2026-09-17
Status: draft, pending user review

## Context

`matpredict detect`'s original search design (see
`docs/superpowers/specs/2026-09-16-mat-detection-pipeline-design.md`,
"Search" section) used `exonerate --model protein2genome` directly
against genomic sequence for two cases: a windowed second pass around
already-found flanking genes, and — after a fix made during
implementation — an unrestricted whole-genome pass for any routed
family with zero hits at all. The whole-genome case is the problem this
revision fixes: `exonerate protein2genome` is a slow, splice-aware
aligner, and running it unrestricted across an entire assembly for
every zero-hit family does not scale.

This revision replaces exonerate's role as a genome-wide search tool
with a cheaper, purpose-built localization step (`tblastn`), and
upgrades the precision step that follows it into a two-tool
cross-validated polish (`miniprot2` and `exonerate --refine region`)
rather than a single tool's unverified output.

## Scope

This spec revises only the **Search** stage of the existing detection
pipeline spec. Everything else (family routing, clustering, scoring,
confidence tiering, idiomorph assignment, the fast path for genomes
with an existing predicted proteome, output format) is unchanged except
where explicitly noted below.

**In scope:**
- Replace the genomic-path search tool from direct `exonerate
  protein2genome` to a `tblastn` localization pass.
- Add a polishing step (`miniprot2` + `exonerate --refine region`) that
  runs only within tblastn-localized candidate regions.
- Add a tool-agreement signal that feeds into existing confidence
  tiering as a bonus, not a gate.
- Export the resulting candidate region (coordinates, per-tool gene
  models, agreement status) as this stage's final output.

**Out of scope (confirmed with the user):**
- Running `augustus`/`helixer` or any other ab initio gene predictor —
  that consumes this stage's exported candidate region and belongs to
  sub-project 3 (gene prediction refinement), not this revision.
- Any change to the fast path (diamond blastp against an existing
  predicted proteome) — that path is untouched; this revision only
  replaces the genomic/fallback path used when no proteome is supplied,
  or when the fast path's second pass would otherwise fall back to
  unrestricted genomic search.
- Retraining or fine-tuning any gene predictor on MAT-specific features
  — the user has flagged that ab initio predictors may perform poorly on
  MAT genes' distinctive features (e.g., extreme brevity for pheromone
  precursors) as a longer-term concern. This is recorded here as
  forward context for sub-project 3/4, not something this spec
  addresses.

## Pipeline stages (revised)

### 1. Localization — tblastn

For a genome with no existing predicted proteome (or whenever the fast
path's per-family second-pass check determines a family's `core_MAT`
gene is missing — see "Search" in the original spec for when the second
pass triggers), run `tblastn` with every routed family's expected
genes as queries — **both `core_MAT` and `flanking_conserved`**,
preserving the dual-anchor principle from the original spec — against
the whole genome. `tblastn`'s translated-nucleotide search has no
splice awareness and gives only approximate hit boundaries, but it is
fast enough to run genome-wide without the performance problem
unrestricted `exonerate protein2genome` had.

`tblastn` hits are converted into the same `SearchHit`-shaped structure
the rest of the pipeline already consumes (family/gene attribution via
the existing record-ID-keyed lookup, not the flat gene-name lookup that
was already fixed as a real bug earlier in this project). These hits
feed into the existing `cluster_hits` gap-based clustering (unchanged)
to produce candidate regions per family — reusing Task 5's clustering
code as-is, just with `tblastn` hits instead of `exonerate` hits as
input.

This stage entirely replaces both of the original spec's genomic-path
uses of `exonerate` (the windowed second pass and the whole-genome
zero-hit fallback) — `exonerate` no longer runs against the whole
genome or an unbounded region under any circumstance after this
revision.

### 2. Polishing — miniprot2 + exonerate --refine region

For each candidate region produced by Stage 1, run both:
- `miniprot2` restricted to the region (padded by a configurable margin
  to avoid boundary truncation), and
- `exonerate --model protein2genome --refine region` restricted to the
  same padded region,

each independently producing a gene model (exon/intron structure,
precise start/end, strand) for the genes expected in that region. Both
tools receive the same curated reference protein(s) as query, and both
run only within the small, tblastn-localized window — never against the
whole genome.

### 3. Agreement as a confidence signal

For each gene, compare `miniprot2`'s and `exonerate --refine`'s
independently-produced boundaries/exon structure:

- **Close agreement** (boundaries and exon/intron structure match
  within a configurable tolerance) becomes an additional positive input
  to the existing confidence tiering (`tiering.assign_tier`, unchanged
  in its core logic) — agreement is a **bonus signal**, not a
  requirement. A family that already reaches "high" tier via the
  original spec's rules (core genes found, flanking evidence where
  applicable) is unaffected by disagreement; agreement can move a
  borderline "medium" case toward "high", consistent with how
  second-pass usage already caps tiers today.
- **Disagreement** does not gate or fail the result. Both tools'
  alignments are kept and reported in the output (see "Output" below);
  the existing tiering logic (core-gene completeness, flanking
  evidence, second-pass usage) still determines the tier on its own
  terms, exactly as it does today — disagreement simply withholds the
  agreement bonus rather than actively penalizing the result.

The exact tolerance/scoring mechanics (how close is "close agreement",
how the bonus is weighted against the existing tier rules) are an
implementation-plan decision, not fixed by this spec — the plan should
propose a concrete, testable rule (e.g., exact coordinate match vs. a
percentage-overlap threshold) and it will be reviewed at that stage.

### 4. Output

This stage's final output for each candidate region is:
- The region's coordinates (contig, start, end — 1-based, fully-closed,
  matching every other coordinate in this project).
- Per-gene evidence from both `miniprot2` and `exonerate --refine`
  (mirroring the existing `GeneEvidence` structure already added to
  `DetectionResult` for the diamond/exonerate fast-path and
  windowed-second-pass evidence — this revision extends that same
  structure to carry both tools' independent gene models rather than
  inventing a new shape).
- An explicit agreement status per gene (agree / disagree / one-tool-only
  — a gene found by only one of the two polishing tools is reported as
  such, not silently dropped).

This output is what sub-project 3 (out of scope here) would consume if
it later runs an ab initio predictor on the same region.

## What this revision removes

- `exonerate`'s unrestricted whole-genome invocation (added as part of
  the "unconditional second pass" fix during the previous implementation
  round) is removed entirely — `tblastn` takes over that role.
- The prior windowed-second-pass mechanism (restricting `exonerate` to a
  region near already-found flanking genes) is superseded by this
  stage's more general localize-then-polish structure — any family
  needing a second look now goes through Stage 1 (localize via
  `tblastn`) → Stage 2 (polish via both tools) uniformly, rather than
  having a separate, narrower flanking-anchored code path.

## Testing / acceptance criteria

- Unit tests for the `tblastn` wrapper, the `miniprot2` wrapper, and the
  revised `exonerate --refine region` invocation all use an injected
  `runner` callable, exactly like every other external-tool wrapper in
  this project — no unit test invokes a live binary.
- A unit test for the agreement-comparison logic using synthetic
  miniprot2/exonerate outputs that agree, and a separate test where they
  disagree, confirming the tiering bonus is applied/withheld correctly
  without affecting the underlying tier rules from the original spec.
- An integration-level test (mirroring the existing
  `test_integration_real_record.py` pattern) confirming the localize
  step's `cluster_hits` reuse works correctly when fed `tblastn`-shaped
  hits instead of `exonerate`-shaped hits.
- Acceptance target: end-to-end run against a genome-only (no supplied
  proteome) test case reaches the same or better confidence outcome as
  the pre-revision unrestricted-exonerate path, but without a
  whole-genome `exonerate protein2genome` invocation ever occurring.

## Out of scope for this spec

- Ab initio gene prediction (augustus/helixer) — sub-project 3.
- Any retraining/fine-tuning of a gene predictor for MAT-gene-specific
  features (e.g., extremely short pheromone-precursor ORFs) — a
  longer-term concern flagged by the user for sub-project 3/4 to
  investigate once ab initio prediction is in scope; recorded here as
  context, not addressed by this revision.
- Any change to the fast path (diamond blastp against a supplied
  proteome) or to family routing, scoring, or idiomorph assignment.
