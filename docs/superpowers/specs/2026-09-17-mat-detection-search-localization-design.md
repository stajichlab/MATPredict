# MATPredict Detection Pipeline: Search Stage Revision — Localize-then-Polish

Date: 2026-09-17 (revised after Fable review)
Status: draft, pending user review

## Context

`matpredict detect`'s original search design (see
`docs/superpowers/specs/2026-09-16-mat-detection-pipeline-design.md`,
"Search" section) used `exonerate --model protein2genome` directly
against genomic sequence for two cases: a windowed second pass around
already-found flanking genes, and — after a bug discovered during
implementation was fixed — an unrestricted whole-genome invocation for
any routed family with zero hits found by other means. The
whole-genome case is the problem this revision fixes: `exonerate
protein2genome` is a slow, splice-aware aligner, and running it
unrestricted across an entire assembly for every zero-hit family does
not scale.

This revision replaces exonerate's role as a genome-wide search tool
with a cheaper, purpose-built localization step (`tblastn`), and
upgrades the precision step that follows it into a two-tool
cross-validated polish (`miniprot` and `exonerate --refine region`,
each run against the same small, tblastn-localized window) rather than
a single tool's unverified output.

**Revision note (second draft):** the first draft was reviewed and
found to have several structural gaps and one factual error, all fixed
below: (1) the tool was misnamed "miniprot2" — the real tool is
`miniprot` (Heng Li); no such thing as miniprot2 exists; (2)
`exonerate --refine region` was described as a region-*restriction*
flag — it is not. `--refine` is an alignment-refinement *strategy*
(re-running exhaustive dynamic programming around a heuristic
alignment); region restriction is done, as it already is elsewhere in
this codebase, by slicing the target FASTA to the padded window first
(`_extract_window`), then running exonerate with `--refine region` on
that slice for improved accuracy within it; (3) the original draft's
"agreement is a tiering bonus" design, its handling of zero/one-tool
outcomes, and its account of what replaces the deleted
`second_pass_used`/Medium-tier mechanism were all underspecified or
internally inconsistent — resolved below with explicit user decisions.

**Revision note (third draft):** implementation surfaced a real gap in
the second draft's Stage 0 design, fixed below: the fast-path rescue was
originally scoped to "zero hits anywhere" (a whole-genome rescue) versus
"a foothold cluster missing one gene" (a narrow, windowed rescue only) —
but this asymmetry was backwards. A family with *partial* annotation
coverage (a receptor gene correctly annotated, its pheromone-precursor
partner missing — this project's own real *M. maydis* mfa1 case) was
getting *less* search than a family with *no* coverage at all, when it
should get the same genome-wide search for its specific missing gene.
**The principle now stated explicitly: a supplied annotation's
completeness for a family's core genes is never trusted, regardless of
how many of that family's genes the annotation did get right** — a
whole-genome annotation missing one real, short, divergent gene while
correctly calling its neighbor is exactly the failure mode `mfa1`
demonstrated (a classic hand-curated locus-specific GenBank deposit
correctly called it; the modern automated whole-genome annotation
dropped it entirely). Stage 0 and Stage 1 below are revised accordingly:
the fast path stays cheap only in the case that actually earns it (every
core gene already found), and any single missing core gene — whether
the family has zero hits or several — triggers the same genome-wide
`tblastn` localization for that specific gene, scoped per (family, gene,
*cluster*) so a real second, independent locus for the same family is
never masked by a different cluster's evidence for the same gene name.

## Scope

This spec revises only the **Search** and the
**second_pass_used**-driven part of the **Boundary calling and
confidence tiering** logic from the existing detection pipeline spec.
Everything else (family routing, clustering's core algorithm, scoring's
fractional-attribution logic, idiomorph assignment, the fast path's use
of diamond against an existing proteome, output format's existing
fields) is unchanged except where explicitly noted below.

**In scope:**
- Replace exonerate's genomic-search role with a `tblastn` localization
  pass for any case that previously ran exonerate against the whole
  genome or an unbounded region.
- Add a polishing step (`miniprot` + `exonerate --refine region`, both
  run against the same sliced, padded window) that runs only within
  tblastn-localized candidate regions.
- Redefine the Medium-tier "second pass" trigger as "gene localized by
  tblastn but not confirmed by either polishing tool" (see "Boundary
  calling" below) — this replaces, rather than adds to, the deleted
  relaxed-exonerate mechanism.
- Report per-gene tool-agreement status (agree / disagree / one-tool-only
  / unpolished) as output data for a human reviewer and for sub-project
  3 — **not** as an input to confidence tiering in v1 (see "Boundary
  calling" below for the reasoning).
- Resolve the fast-path/genome-only interaction explicitly (see
  "Pipeline stages", stage 0).

**Out of scope (confirmed with the user):**
- Running `augustus`/`helixer`/`braker` or any other ab initio gene
  predictor — that consumes this stage's exported candidate region and
  belongs to sub-project 3 (gene prediction refinement), not this
  revision.
- Feeding tool-agreement into confidence tiering in v1 (see "Boundary
  calling" below) — deferred until the leave-one-out benchmark can show
  agreement actually correlates with correctness, since both polishing
  tools receive the same reference protein and the same window and are
  not truly independent evidence sources (they will agree wherever
  local alignment is unambiguous and diverge only where homology is
  genuinely weak — this measures alignment clarity for one gene model,
  not locus-level correctness).
- Retraining or fine-tuning any gene predictor on MAT-specific features
  — the user has flagged that ab initio predictors may perform poorly on
  MAT genes' distinctive features (e.g., extreme brevity for pheromone
  precursors) as a longer-term concern for sub-project 3/4;  recorded
  here as forward context, not addressed by this revision.

## Pipeline stages (revised)

### Stage 0 — routing between the fast path and genome-wide localization

This revision applies in two different situations that must be handled
differently, not identically (a gap in the first draft), and the
fast-path case has since been narrowed to reflect the third-draft
principle above (an annotation's completeness for a family's core genes
is never trusted):

1. **Genome-only path** (no predicted proteome supplied, or the
   exhaustive no-taxid case): no prior search has happened yet, so
   there is no existing cluster to anchor a window to. This goes
   through the full Stage 1 (genome-wide `tblastn` localization) before
   Stage 2 (polish).
2. **Fast path, no rescue needed** (a proteome was supplied and diamond
   found every one of a family's `core_MAT` genes): this is the ONLY
   case that stays cheap and skips Stage 1 entirely — the fast path's
   performance benefit is real, but it is earned by the annotation
   actually being complete for that family, not assumed.
3. **Fast path, rescue needed** (a proteome was supplied, but diamond
   left one or more of a family's `core_MAT` genes missing — whether
   the family has an existing foothold cluster with one gene still
   missing, or has zero hits for the family at all): both cases now
   receive the SAME treatment — the specific missing gene(s)' curated
   proteins are added to the same batched, genome-wide `tblastn`
   localization call the genome-only path already runs (one call per
   `run_pipeline` invocation, not one per family or per gene). This
   replaces the second draft's asymmetric design, where a zero-hit
   family got genome-wide search but a partial-foothold family only got
   a narrow window — that asymmetry produced *less* search for *more*
   annotation coverage, backwards from the intent.

**Cluster-aware rescue eligibility.** Whether a specific gene is
"missing and needs rescue" is evaluated per `(family, gene, cluster)`,
not per `(family, gene)` alone. A family can have more than one real,
independent locus in one genome (tetrapolar species with unlinked loci;
homothallic/heterothallic switching-cassette species) — if gene X is
found in cluster A but missing from cluster B, cluster B's copy of gene
X still needs its own rescue, even though the family as a whole "has"
gene X somewhere. Scoping eligibility to the specific cluster, not just
the family, is what makes this catch a real second-locus gap instead of
treating the family-wide presence of a gene as evidence that every
cluster of that family has it.

The genome-wide `tblastn` call itself costs nothing extra either way:
it always queries the family's full expected gene set (`core_MAT` and
`flanking_conserved`) against the whole genome in one invocation,
regardless of which specific genes/clusters are eligible for rescue —
`families`/eligibility only controls which of the resulting hits are
kept after the fact. Widening eligibility from per-family to
per-cluster is therefore a question of how strict that post-hoc
filtering should be (accepting a modest increase in spurious-hit
surface, mitigated by the existing scoring/clustering/ambiguity-floor
thresholds a spurious hit must still clear before being promoted to a
called locus), not a batching or performance tradeoff.

Clustering runs **once** for the genome-only path: `tblastn` hits for
every routed family (both `core_MAT` and `flanking_conserved` genes,
preserving the dual-anchor principle) are batched into a single hit
pool alongside any other search-stage hits, and `cluster_hits`
(unchanged) runs once over that pool — mirroring the existing
batched-exonerate-call pattern for the zero-hit case, now applied to
`tblastn`. For the fast-path rescue case, a rescued hit that lands
within an existing cluster's span is folded into that cluster directly
(without re-running `cluster_hits`, the same mechanism the original
windowed second pass used); a rescued hit that lands elsewhere becomes
the seed of its own new cluster, since it may represent a genuine
second, independent locus for that family.

### Stage 1 — Localization (genome-only path only)

Run `tblastn` once, batched across every routed family's expected
genes (both `core_MAT` and `flanking_conserved`) as queries, against
the whole genome (one `makeblastdb`-equivalent database build plus one
`tblastn` invocation, not one per family). Disable the default
low-complexity filter (`-seg no`) — the default SEG filter would
suppress exactly the short, simple pheromone-precursor queries this
project's own domain history (sub-project 1's *M. maydis* mfa1 case)
already identified as a detection blind spot.

`tblastn` HSPs are normalized before becoming `SearchHit`s:
minus-strand HSPs report `sstart > send` and must be normalized to
`start <= end` with `strand` derived from the HSP's frame sign — stated
explicitly here since this is a real, previously-hit bug class in this
project (coordinate/offset handling for a new search tool). Multiple
HSPs from introns in one protein-vs-genome alignment are expected and
are correctly merged by `cluster_hits`'s gap-based grouping; no dedupe
is needed at this stage since `score_cluster` counts distinct gene
*names*, not distinct HSPs.

This stage entirely replaces the original spec's unrestricted
whole-genome `exonerate` fallback, and now runs for BOTH the
genome-only path and any fast-path family needing rescue (per Stage 0's
third-draft revision) rather than only for zero-hit families. The
original spec's flanking-anchored relaxed search is superseded by
Stage 2's real polish-with-both-tools step, which now runs uniformly on
every localized/rescued cluster rather than a separate "relaxed
exonerate" mode — see "Boundary calling" for what replaces the
confidence-tier effect that mechanism used to have.

### Stage 2 — Polishing (both paths)

For each candidate window — either a Stage 1 cluster's span (genome-only
path) or an existing cluster's span (fast-path rescue), both padded by
a configurable margin — run **both**:

- `exonerate --model protein2genome --refine region`, and
- `miniprot`,

against the same sliced window (via `_extract_window`, the existing
FASTA-slicing helper, with the same coordinate re-basing already used
for the windowed second pass — reused unchanged, since padding
introduces the same offset-tracking requirement that mechanism already
solved). Both tools receive the same curated reference protein(s) as
query.

The padding margin's default is derived from each candidate query
protein's own curated length (e.g., a multiple of the protein's
nucleotide-equivalent length, plus a max-intron allowance) rather than
a single fixed constant — a fixed default risks truncating a real gene
when `tblastn` only hit a short conserved domain (common for divergent
MAT proteins, e.g. an HMG or homeodomain). The exact multiplier is an
implementation-plan decision but must be named and configurable, not
silently hardcoded.

When a family/gene has multiple curated candidate proteins, run polish
against each and select "the" model per gene per tool by a named,
tool-appropriate score (bitscore or aligned length — **not** raw
percent identity, since `tblastn`'s `pident`, `miniprot`'s identity,
and `exonerate`'s identity are not directly comparable numbers across
tools) before comparing tool-vs-tool for agreement.

### Stage 3 — Per-gene status and canonical coordinates (replaces "agreement as confidence signal")

For each gene in a candidate window, after Stage 2, classify it into
exactly one status:

- **`polished_agree`**: both `exonerate --refine` and `miniprot`
  produced a gene model, and their boundaries/exon structure agree
  within a configurable tolerance.
- **`polished_disagree`**: both tools produced a model, but they
  disagree beyond tolerance. Both models are kept in the output (never
  silently dropped) — this is informational for a human reviewer (real
  exon-structure disagreement between two independent aligners is also
  a known indicator of pseudogenization/frameshift, which is directly
  relevant to MAT loci) but does **not** affect confidence tiering.
- **`polished_single`**: exactly one of the two tools produced a usable
  model (the region may be too short, degraded, or genuinely divergent
  for the other tool). That tool's model is used for canonical
  coordinates.
- **`unpolished`**: neither polishing tool produced a usable model, but
  `tblastn` did localize a rough hit. This status is what feeds
  confidence tiering below — it **replaces** the deleted
  `second_pass_used` mechanism entirely, keeping the same effect the
  original spec intended (a gene confirmed only by the less precise
  search method caps the tier at Medium, never High) without needing a
  "standard vs relaxed" search-parameter distinction, since there is no
  such distinction left in this design.

**Canonical coordinates for output/scoring**: `polished_agree` and
`polished_disagree` both use `exonerate --refine`'s coordinates as
canonical (it is the more directly integrated tool, already producing
protein-based, intron-aware alignments matching sub-project 1's own
curation conventions) — `miniprot`'s model is retained as the
corroborating/comparison evidence in the output, never discarded.
`polished_single` uses whichever tool succeeded. `unpolished` uses the
raw `tblastn` HSP's coordinates (already normalized per Stage 1),
explicitly flagged as approximate in the output.

**Effect on `score_cluster`/tiering**: a gene in any of the four
statuses above counts as "found" for `score_cluster`'s fractional
scoring (an `unpolished` gene is still real evidence the gene exists,
just imprecisely located) — the distinction only affects the
confidence tier, not whether the gene counts toward `fraction_found`.
`assign_tier`'s existing High/Medium/Low structure is otherwise
unchanged: `unpolished` for any of a family's genes caps that family's
tier at Medium, exactly where `second_pass_used` did before; a family
with every gene `polished_agree`/`polished_disagree`/`polished_single`
reaches High under the same core/flanking-evidence rules the original
spec already defined. **Tool agreement itself is not consulted by
`assign_tier`** — `polished_agree` and `polished_disagree` are treated
identically for tiering purposes (both are "polished"), per the
decision that agreement is not yet validated as a correctness signal.

### 4. Output

`GeneEvidence` (already added to `DetectionResult` for the
diamond/exonerate fast path) is extended — this is a genuine shape
change, not a same-shape reuse, and the implementation plan must treat
it as such — to carry, per gene:
- the canonical coordinates and the tool that produced them,
- the per-gene status (`polished_agree` / `polished_disagree` /
  `polished_single` / `unpolished`),
- both tools' individual models when both ran (even if only one is
  canonical), so a `polished_disagree` case is fully inspectable.

The candidate window's own coordinates (contig, start, end — 1-based,
fully-closed, matching every other coordinate in this project) are
reported alongside the per-gene evidence; whether the window's reported
span is the raw `tblastn` cluster span or extends to cover any
polished model that reached past it (via padding) must be stated
explicitly by the plan — the polished model's own coordinates are
authoritative for the gene; the window's reported span should cover the
union of all its genes' canonical coordinates, not be silently
truncated to the original `tblastn` HSP extent.

This output is what sub-project 3 (out of scope here) would consume if
it later runs an ab initio predictor on the same region.

## What this revision removes

- `exonerate`'s unrestricted whole-genome invocation (added as part of
  the "unconditional second pass" fix during the previous implementation
  round) is removed entirely — `tblastn` takes over that role for the
  genome-only path.
- The prior windowed-second-pass mechanism's "relaxed exonerate
  parameters" (`--percent`/`--score` threshold changes) are removed —
  every fast-path rescue case now runs the same real
  polish-with-both-tools step described above, not a relaxed
  single-tool retry.
- The second draft's asymmetric rescue scope (genome-wide search for a
  zero-hit family, but only a narrow window for a partial-foothold
  family) is removed — both now receive the same genome-wide `tblastn`
  rescue for their specific missing gene(s), scoped per
  `(family, gene, cluster)` per the third-draft revision above.
- `second_pass_used` as a tiering input is removed and replaced by the
  `unpolished` status (Stage 3) — same tiering *effect* (caps at
  Medium), different, more precisely-defined trigger.
- `SearchHit.method` values `exonerate_genome`/`exonerate_genome_relaxed`
  are retired; new values are needed for `tblastn_genome`,
  `miniprot_genome`, and `exonerate_refine`. `report.py`'s GFF3/YAML
  output already emits `method` into user-facing output — the
  implementation plan must audit every place a `method` string is
  matched, displayed, or tested, not just the `search.py` producer side.
- `NotDetectedFamily`'s literal reason text referencing "the relaxed
  whole-genome second pass" must be updated to describe the new
  localize-then-polish flow instead.
- The original spec's "Fallback path" paragraph (naming "exonerate or
  minimap2") is superseded by this document for the genome-only case.

## Testing / acceptance criteria

- Unit tests, all using injected `runner` callables per this project's
  existing convention (no unit test invokes a live binary):
  - `tblastn` minus-strand HSP normalization (start/end swap, strand
    derived from frame).
  - Window-offset re-basing for **both** polishing tools (a previously
    real bug class in this project — exonerate's windowed second pass
    already needed this fix once).
  - All four Stage 3 statuses (`polished_agree`, `polished_disagree`,
    `polished_single`, `unpolished`), including that `unpolished`
    caps a family at Medium tier and the other three do not, and that
    `polished_agree` vs `polished_disagree` produce identical tiering
    outcomes (proving agreement is genuinely not consulted by
    `assign_tier`, not just untested).
  - The best-model-per-gene-per-tool selection logic (multiple curated
    candidate proteins, correct score used, not raw cross-tool
    identity).
  - The fast path's cheap case (Stage 0 case 2): confirm a family with
    every core gene found by diamond skips Stage 1 entirely and never
    triggers `tblastn`.
  - The fast path's rescue case (Stage 0 case 3), both shapes: confirm
    a zero-hit family AND a partial-foothold family both trigger the
    same batched genome-wide `tblastn` call for their specific missing
    gene(s) — not a narrower, windowed-only search for the
    partial-foothold case.
  - Cluster-aware rescue eligibility: confirm that when the same family
    has gene X found in one cluster but missing from a second,
    independent cluster, the second cluster's copy of gene X still gets
    rescued — eligibility must be keyed per `(family, gene, cluster)`,
    not per `(family, gene)` alone (a bug found and fixed during this
    revision's final review: family-wide presence of a gene was
    incorrectly treated as covering every cluster of that family).
  - The genome-only/rescue single-batched-`tblastn`-call behavior
    (confirm one call covers every routed family's genes needing
    localization, not one call per family or per gene) — assert
    directly on the injected runner's captured invocation that
    `--target`/equivalent is never the full `genome_fasta` for the
    *polishing* tools (only the sliced window is ever passed to them),
    and that exactly one `tblastn` invocation covers the whole run.
  - Multi-gene window handling: both polishing tools must correctly
    handle a window containing more than one gene's alignment (the
    normal case for a real MAT locus, not an edge case) — group each
    tool's output by its own alignment record, filter to the
    specifically requested gene, and select the best surviving record
    by the tool's own score column, never by comparing identity across
    tools (a bug found and fixed during this revision's final review:
    both wrappers originally kept only the first feature record in the
    tool's output, so a multi-gene window let one gene's data silently
    contaminate or crowd out another's).
- Integration-level test (mirroring
  `tests/detect/test_integration_real_record.py`) confirming the
  genome-only localize-then-polish flow works end-to-end against a real
  curated record with a real, on-disk sliced-window FASTA (not just
  mocked coordinates), covering the offset re-basing for real.
- Acceptance targets:
  - Run the existing leave-one-out benchmark (sub-project 6) against
    `assembly`-type curated records (the ones with a real genome, not
    just a locus fragment) for gene recall and boundary accuracy,
    comparing this revision's results to a recorded pre-revision
    baseline run — not a single-tier pass/fail check, since a
    three-value tier can't distinguish a real regression from a lucky
    pass.
  - A measured wall-clock bound on a named real test assembly,
    recorded and compared against the pre-revision unrestricted-exonerate
    baseline on the same assembly — since eliminating that runtime cost
    is this revision's entire motivation, the acceptance criteria must
    include a number, not just "faster in principle."

## Out of scope for this spec

- Ab initio gene prediction (augustus/helixer/braker) — sub-project 3.
- Feeding tool-agreement into confidence tiering — deferred pending
  leave-one-out benchmark evidence that it correlates with correctness.
- Any retraining/fine-tuning of a gene predictor for MAT-gene-specific
  features (e.g., extremely short pheromone-precursor ORFs) — a
  longer-term concern flagged by the user for sub-project 3/4 to
  investigate once ab initio prediction is in scope; recorded here as
  context, not addressed by this revision.
- Any change to family routing, `score_cluster`'s fractional-attribution
  algorithm, or idiomorph assignment.

**Forward note for sub-project 3 (not a requirement of this spec):**
the user has indicated sub-project 3 should be able to run a *targeted*
ab initio predictor (augustus/helixer/braker) restricted to a candidate
region this stage identifies, rather than genome-wide. This spec's
Stage 4 output (a candidate region's coordinates plus its per-gene
evidence) is the natural input for that — no change is needed here to
support it, but sub-project 3's own design should confirm this output
shape is sufficient before adding anything new to it.
