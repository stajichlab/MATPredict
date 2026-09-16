# MATPredict Sub-project 2: Detection/Annotation Pipeline

Date: 2026-09-16 (revised after Fable review)
Status: draft, pending user review

## Context

Sub-project 1 (the curated MAT locus reference database, spec at
`docs/superpowers/specs/2026-09-16-mat-reference-database-design.md`)
has produced 30 accepted records across Mucoromycota, Basidiomycota,
and Ascomycota, spanning roughly 10 distinct locus architectures
(bipolar 2-flanking-gene, bipolar-fused, tetrapolar HD/PR, 4-sublocus
tetrapolar, fused homothallic PM, switchable mat1/mat2/mat3, MATa1/
alpha1/alpha2, HMG-domain MATA/MATB, tetrapolar a/b smut loci). This
sub-project builds the tool that consumes that reference data to find
and characterize MAT loci in a new, uncurated genome.

Per the original project scope: "given a genome, identify the MAT
locus, with or without annotation already provided on the genome, and
with or without classification of what taxonomy id it has." Background
knowledge established at project start (and reconfirmed during
sub-project 1 curation — the *Mycosarcoma maydis* mfa1 case, a 41-aa
pheromone-precursor ORF missing from automated whole-genome
annotation): general-purpose gene predictors often do poorly on the
core MAT genes themselves, while flanking genes are frequently easier
to detect. This sub-project's search strategy is built around that
asymmetry, **but see the explicit limitation in "Short pheromone-
precursor genes" below**: the standard homology methods this spec
proposes cannot actually detect genes as short as mfa1 either, and the
spec says so rather than implying otherwise.

**Revision note**: this is the second draft. The first draft was
reviewed and found to have several structural gaps, checked directly
against the real accepted records rather than the reference DB spec's
abstractions: confidence tiers that most families could never reach,
underspecified family routing given how many distinct `locus_name`
entries now exist per phylum, no handling for cross-family homology
contamination or multi-locus genomes, a boundary-window heuristic that
doesn't hold across the DB's actual size range (809 bp to 148 kb), a
leave-one-out benchmark design that would leak or trivially pass for
most of the current DB, and no idiomorph-assignment stage despite
classification being a core project goal. All of these are addressed
below.

## Scope

Build `matpredict detect`, a new CLI subcommand that:
1. Accepts a genome (FASTA, optionally with existing gene predictions)
   and an optional NCBI taxid.
2. Selects candidate curated reference families using taxonomy
   (narrowed as far as the taxid's lineage and each `order.yml` entry's
   declared taxonomic scope allow — not just phylum-level).
3. Searches for MAT locus genes via homology, using a **dual-anchor**
   strategy (`core_MAT` and `flanking_conserved` genes are both
   independently search targets).
4. Clusters hits per contig/region, calls locus boundaries, assigns a
   confidence tier per family definition (not one fixed rule that only
   some families can reach), and reports **all** plausible loci found
   (a genome can genuinely contain more than one MAT-related locus).
5. Assigns idiomorph where the matched family's vocabulary supports it.
6. Emits a GFF3 of each predicted locus plus a validation-style report.

This sub-project also designs (per project scope, "alongside 1 and 2")
the Sn/Sp benchmark suite (sub-project 6), implemented as a
species/genus-level leave-one-out test against sub-project 1's own
curated records, with explicit sample-size reporting per family.

**Explicitly out of scope for this sub-project**: profile HMM
construction (deferred until a `locus_name` family has enough curated
diversity to support one), *de novo* gene structure prediction via
BRAKER/AUGUSTUS/Helixer (sub-project 3), ML/LLM-based sequence
refinement (sub-project 4), RNA-seq integration, and resolving which
copy of a switching-system's cassettes (mat1 vs mat2/mat3, MAT vs
HML/HMR) is transcriptionally active (homology alone cannot tell —
see "Multiple loci and switching systems" below).

## Why pairwise homology search, not profile HMMs, for v1

Most `locus_name` families in the current reference DB have only 1-3
curated examples (e.g. `MATyl` has exactly 1; `Balpha`/`Bbeta` have 1
each). A profile HMM built from a single sequence, or from 2-3 closely
related strains of the same species, would not capture real
across-species sequence diversity and would effectively just be an
expensive way to do a single-sequence search — while risking false
confidence from a formal-looking model. Pairwise search (diamond/blastp
for protein-vs-protein, exonerate/minimap2 for protein-vs-genome)
degrades gracefully to exactly one comparison per curated example and
requires no minimum family size. As families grow (particularly `MAT`,
which already has the most examples across Pezizomycotina/
Mucoromycota/Basidiomycota — noting `MAT` as a `locus_name` string is
reused across all three phyla with completely different gene sets in
each; **the real family key is `(phylum, locus_name)`, never
`locus_name` alone**), building real profile HMMs for those specific
families becomes a natural, separate follow-up — not a v1 blocker.

## Family routing: taxonomy narrows candidates, doesn't always pick one

`order.yml` needs a new field per locus entry to make routing concrete:
`taxonomic_scope`, an explicit list of taxids (or a single taxid
representing the scope's common ancestor) that entry applies to. For
example, Ascomycota's `order.yml` currently has 8 `locus_name` entries
(`MAT`, `PM`, `mat1`, `mat2`, `mat3`, `MATsc`, `MATyl`, `MTL`) with
disjoint real taxonomic ranges (`MATsc` → Saccharomycetaceae; `mat1` →
Schizosaccharomycetaceae; `PM` → Pneumocystidales/Taphrinales; `MAT` →
Pezizomycotina) — phylum-level narrowing alone still leaves 8
candidates to try for any Ascomycota genome. `taxonomic_scope` lets
routing narrow to the entries whose scope actually contains (or is an
ancestor/descendant of) the input taxid's lineage.

**When a taxid is given**: resolve its lineage via taxonkit (reusing
`MATPredict.db.taxonomy.resolve_lineage`), then select every
`(phylum, locus_name)` entry across all phyla whose `taxonomic_scope`
overlaps that lineage. This is normally one phylum's entries, narrowed
further by scope — but not always exactly one entry (see "Multiple
loci and switching systems" below for tetrapolar/multi-cassette
species that legitimately need several).

**When no taxid is given, or its lineage matches no declared scope**:
fall back to attempting every `(phylum, locus_name)` family in the
database. This is the exhaustive path — see "Cross-family
contamination and ambiguous matches" for how results from this path
are scored and reported (never a silent single winner).

`db/_schema/order.schema.yaml` gains a required `taxonomic_scope` field
per locus entry; populating it for the entries built during sub-project
1 (using each accepted record's own resolved taxid/lineage as the
scope) is an implementation task for this sub-project, not sub-project
1 — sub-project 1's `order.yml` files predate this requirement.

## Pipeline stages

### 1. Search

**Fast path (existing gene predictions available)**: search the
genome's predicted proteins directly against the curated reference
protein set for each candidate family (diamond or blastp). **Any
expected `core_MAT` gene not found among the proteome hits within a
plausible window always triggers the genomic second-pass search below
— this check is unconditional, not gated on whether flanking genes
were found first.** This closes the exact blind spot sub-project 1
discovered in existing genome annotations (mfa1-style short/divergent
genes silently absent from a provided proteome) — the whole point of
this pipeline is to not simply trust that an existing annotation is
complete.

**Fallback path (genome FASTA only, no predictions), and the
genomic second pass above**: use spliced protein-to-genome alignment —
exonerate (`--model protein2genome`) or minimap2 (`-x splice`) — to
align each curated reference protein directly against the genomic
sequence. Naive six-frame translation was considered and rejected:
curated records already show many intron-containing MAT genes (e.g.
Cryptococcus SXI1, Coprinopsis HD1), which six-frame translation would
fragment or miss entirely.

**Dual-anchor targets**: both `core_MAT` and `flanking_conserved` genes
from each candidate family are searched independently in the same
pass — neither is a prerequisite for searching the other, and **many
families (all of Basidiomycota's currently-accepted records, PM,
mat1/mat2, MATsc, bLocus) have zero `flanking_conserved` genes in
their curated gene list at all** — the confidence-tiering logic below
does not assume every family has flanking genes to anchor on.

### 2. Clustering, cross-family scoring, and ambiguity

All hits (both anchor types, across all attempted families) are first
grouped into spatial clusters: hits within a configurable maximum
intergenic gap (default ~25 kb, tunable — not derived from any single
family's locus span, since spans range from ~800 bp to ~148 kb even
*within* one family) on the same contig are one cluster, regardless of
strand (real curated loci mix strands within one locus — e.g. *Y.
lipolytica*'s MATA has genes on both strands — so strand consistency is
never a hard clustering requirement, only an input to orientation
scoring below).

Each cluster is then scored **per candidate family independently**, as
the fraction of that family's expected genes (per its `order.yml`
entry) found in the cluster, not a summed bit-score — summed score
would systematically favor families with more genes regardless of
whether they're the right match (e.g. an 8-gene family scoring higher
than a correct 2-gene family purely on gene count). Gene-order/
orientation agreement with the matched record's `reference_orientation`
is a secondary tie-breaker, not the primary score.

**Cross-family contamination is a structural risk, not an edge case**:
homeodomain genes (HD1/HD2, bE/bW, Z/Y, SXI1/SXI2, MATa1/alpha2), HMG
genes (MAT1-2-1, mata2), and pheromone receptors (STE3, pra1, bar3,
bbr2) recur across many families with real sequence similarity to each
other. When multiple families score above a floor threshold for the
same cluster, **report all of them as candidates with an explicit
"ambiguous" flag** — never silently pick a single winner by score
alone. This applies in both the taxonomy-narrowed and exhaustive
search paths.

### 3. Boundary calling and confidence tiering

Tiering is defined **per family based on what evidence that family can
actually produce**, not one fixed rule:

- **High**: for a family with `flanking_conserved` genes in its
  `order.yml` entry — all expected `core_MAT` genes found plus at
  least one `flanking_conserved` hit in the same cluster. For a family
  with *no* `flanking_conserved` genes at all (a large fraction of the
  current DB) — all expected `core_MAT` genes found in the cluster,
  in an order/orientation consistent with `reference_orientation`
  (multi-gene internal consistency substitutes for flanking evidence
  when no flanking genes exist to check against).
- **Medium**: some but not all expected genes found; or a `core_MAT`
  gene found only via the relaxed genomic second-pass (see below); or
  (for flanking-bearing families) flanking genes found and clustered
  but the corresponding `core_MAT` gene remains unconfirmed even after
  the second pass.
- **Low**: a single gene hit with no other expected genes from the
  same family found nearby, or genuinely no cluster clears the
  minimum floor for any attempted family — reported as **"not
  detected"**, explicitly listing which families were attempted and
  why each fell short (never silently omitted).

**The flanking-anchored second pass, restated concretely**: when
`flanking_conserved` hits cluster together but the expected `core_MAT`
gene from the same family doesn't hit at the standard threshold, run a
second, more sensitive search restricted to the genomic window
between/around the flanking hits (relaxed e-value threshold and/or
exonerate's more permissive alignment modes) before concluding the
core gene is absent.

### 4. Short pheromone-precursor genes: a stated limitation, not a silent gap

Genes like *U. maydis* mfa1 (41 aa) or the various `bbp`/`bap`
pheromone precursors in this database (60-90 bp CDS) are **not
reliably detectable by diamond, blastp, or exonerate at standard
significance thresholds** — sequences this short don't produce
statistically distinguishable alignment scores from noise with
general-purpose homology tools. This is the pipeline's own motivating
example from sub-project 1, and this spec does not claim to solve it
with the methods above. Two options, not mutually exclusive:
(a) explicitly report these specific genes as **"not searchable by
this method"** in the output rather than silently reporting them
absent, so a human reviewer knows the gap is a tool limitation, not a
negative result; (b) add a narrowly-scoped short-ORF/motif rescue scan
(restricted to an already-called locus window, only run after
boundaries are established from other genes) that looks for ORFs in
the expected size range with composition similar to the curated
examples, flagged as low-confidence candidates requiring manual
review, never auto-confirmed. Option (b) is a small, explicit v1
component given how central this problem is to the project's original
motivation — not deferred to sub-project 4 wholesale, though any
deeper motif/ML-based refinement of it is.

### 5. Multiple loci and switching systems

A genome can genuinely contain more than one MAT-related locus:
tetrapolar species have 2-4 unlinked loci (the DB already has HD+PR
for *Coprinopsis*, `Aalpha`/`Balpha`/`Bbeta` for *Schizophyllum*,
`aLocus`/`bLocus` for *Ustilago*), and switching-system species carry
homologous silent cassettes alongside the active locus (*S. pombe*'s
mat1/mat2/mat3, *S. cerevisiae*'s MAT/HML/HMR). The pipeline reports
**one result per spatial cluster**, not one result per genome — when
taxonomy narrows to a multi-locus family set (e.g. both `HD` and `PR`
are in scope for the input taxon), the pipeline actively searches for
all of that clade's expected loci and reports them together as a set,
not just the first one found.

**Switching-system cassettes cannot be distinguished as active vs.
silent by homology search alone** — a silent donor cassette is, by
definition, sequence-similar to the active locus. When multiple
clusters match the same switching family (`mat1`/`mat2`/`mat3` or
`MATsc`+its HML/HMR-equivalent, if curated), report all of them with
expression status **"undetermined"** rather than guessing which is
active. Resolving this would need synteny context or expression data,
both out of scope here.

### 6. Fragmented assemblies

When a family's expected genes are found on different contigs (a real
possibility with a fragmented assembly), the pipeline reports a
multi-segment locus using the same `segments[]` structure sub-project
1's schema already defines for this purpose, with `contig_edge_distance`
populated where a hit sits near a contig boundary. A multi-segment call
is downgraded one confidence tier from what it would otherwise score,
reflecting the added uncertainty about the true intervening sequence.

### 7. Idiomorph assignment

For enum-vocabulary families (`order.yml`'s `vocabulary_type: enum`,
e.g. `MAT1-1`/`MAT1-2`, `Plus`/`Minus`, `alpha`/`a`): idiomorph is
assigned from which idiomorph-restricted genes were actually found
(matching `order.yml`'s `present_in_idiomorphs` mapping for each gene).
For pattern-vocabulary families (`HD`/`PR`-style multiallelic loci,
e.g. `A43`, `B43`): **homology search cannot call the specific allele
number** — output reports the locus type/family as confirmed but the
specific allele as **"undetermined"**, never guessed from which
curated example happened to score best.

### 8. Output

- **GFF3**: predicted locus region(s) and gene features, using the
  same 1-based fully-closed coordinate convention and `role`/`present`
  semantics as sub-project 1's schema, so predictions and curated
  records are directly comparable. Multi-segment loci use multiple
  `##sequence-region` blocks per the schema's existing `segments[]`
  model.
- **Validation-style report**: per detected cluster — matched
  family/families (with the ambiguous flag when more than one scores
  above floor), per-gene identity/coverage against the matched
  reference protein(s), which specific curated record(s) the match was
  scored against, confidence tier, idiomorph call (or "undetermined"),
  and which expected genes were not found, explicitly distinguishing
  "not found" from "not searchable by this method" (short-ORF case).
  Deliberately mirrors `metadata.yaml`'s `validation` block shape from
  sub-project 1.

## Sn/Sp benchmark suite (sub-project 6)

**Method**: species/genus-level leave-one-out against sub-project 1's
own curated records — not per-record leave-one-out. Per-record
holdout leaks: several families hold Plus/Minus (or similarly paired)
idiomorph records from the *same species*, so leaving out one still
leaves near-identical flanking sequence from its pair in the reference
set, trivially inflating recovery. Holding out at species (or genus,
for very close relatives) level removes that leakage.

**What's actually being tested**: most current records
(`coordinate_provenance: insdc_nucleotide`) are locus-only sequence
fragments (1-20 kb), not full genome assemblies — running detection
against the fragment itself can measure whether the pipeline recovers
the correct genes from a held-out species, but **cannot** test false-
positive rate or boundary over-extension, since there's no surrounding
genomic context to falsely match. Where a full assembly exists for a
held-out record (the `assembly`-type `sequence_source` records), run
detection against the full assembly instead of the fragment specifically
to get a real false-positive/fragmentation test; where only a fragment
exists, report gene-identity recall only and say so in the output.

**Reporting**: Sn/Sp computed per `(phylum, locus_name)` family, with
the **post-holdout reference-set size printed alongside every number**.
Families with 2 or fewer total examples produce "n/a — insufficient
data" after holdout removes one, rather than a misleadingly precise
score from an n=0 or n=1 reference set. The database currently has
zero `present: false` (curated absence) genes, so specificity for gene
*absence* calls is not yet testable — flagged as a known data gap for
ongoing sub-project 1 curation to address (a curated example with a
confirmed-absent expected gene would be genuinely useful test data),
not something this benchmark can fabricate.

Records marked `excluded_from_coordinate_benchmark: true` (the *Z.
rouxii* translocation-background record, the *Taphrina* tblastn-only
record) are excluded from boundary-accuracy scoring specifically, but
remain usable for gene-identity-only scoring.

## Testing / acceptance criteria for this sub-project

- Unit tests for taxonomy/scope-based routing (including the
  multi-candidate-family case), clustering (max-gap, cross-strand),
  per-family scoring and the ambiguous-match flag, all three confidence
  tiers (including flank-less-family "high" and the relaxed second-pass
  path), multi-segment/fragmented-assembly handling, and idiomorph
  assignment (both enum and "undetermined"-pattern cases) — using
  mocked search-tool output, not live diamond/exonerate runs, matching
  sub-project 1's established testing pattern (never call live external
  tools from unit tests).
- The leave-one-out benchmark itself, run against the real curated DB,
  as an integration-level check (not a fast unit test) — reports Sn/Sp
  per family with sample size shown, not a pass/fail gate, since most
  families are expected to report "n/a" at the current DB size until
  more curation happens.
- Acceptance target: the pipeline runs end-to-end (routing → search →
  clustering/tiering → idiomorph assignment → GFF3 + report output)
  against at least one real held-out record from each of the three
  phyla, correctly reports "ambiguous" for at least one deliberately
  constructed cross-family test case, and the leave-one-out benchmark
  produces a real, non-fabricated Sn/Sp report (including its "n/a"
  entries) covering every family in the database.

## Out of scope for this sub-project

- Profile HMM construction for any `locus_name` family
- *De novo* gene structure prediction (BRAKER/AUGUSTUS/Helixer) --
  sub-project 3
- ML/LLM-based boundary/motif refinement beyond the narrow short-ORF
  rescue scan described above -- sub-project 4
- RNA-seq integration
- Resolving active-vs-silent cassette status in switching systems
- A public web resource for browsing detection results -- sub-project 5
