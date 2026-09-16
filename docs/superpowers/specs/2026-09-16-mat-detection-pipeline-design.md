# MATPredict Sub-project 2: Detection/Annotation Pipeline

Date: 2026-09-16
Status: draft, pending user review

## Context

Sub-project 1 (the curated MAT locus reference database, spec at
`docs/superpowers/specs/2026-09-16-mat-reference-database-design.md`)
has produced 30 accepted records across Mucoromycota, Basidiomycota,
and Ascomycota, spanning roughly 10 distinct locus architectures
(bipolar 2-flanking-gene, bipolar-fused, tetrapolar HD/PR, 4-sublocus
tetrapolar, fused homothallic PM, switchable mat1/mat2/mat3, MATa1/
alpha1/alpha2, HMG-domain MATA/MATB). This sub-project builds the tool
that consumes that reference data to find and characterize MAT loci in
a new, uncurated genome.

Per the original project scope: "given a genome, identify the MAT
locus, with or without annotation already provided on the genome, and
with or without classification of what taxonomy id it has." Background
knowledge established at project start (and reconfirmed during
sub-project 1 curation, see the *Mycosarcoma maydis* mfa1 case — a
41-aa pheromone-precursor ORF missing from automated whole-genome
annotation): general-purpose gene predictors often do poorly on the
core MAT genes themselves, while flanking genes are frequently easier
to detect. This sub-project's search strategy is built around that
asymmetry.

## Scope

Build `matpredict detect`, a new CLI subcommand that:
1. Accepts a genome (FASTA, optionally with existing gene predictions)
   and an optional NCBI taxid.
2. Selects the appropriate curated reference protein set (narrowed by
   taxonomy when available, exhaustive search when not).
3. Searches for MAT locus genes via homology, using a **dual-anchor**
   strategy (core_MAT genes and flanking_conserved genes are both
   independently search targets, not a fixed core-first search order).
4. Calls locus boundaries from clustered hits and reports a
   three-tier confidence level.
5. Emits a GFF3 of the predicted locus plus a validation-style report.

This sub-project also designs (per project scope, "alongside 1 and 2")
the Sn/Sp benchmark suite (sub-project 6), implemented as a leave-one-out
test against sub-project 1's own curated records.

**Explicitly out of scope for this sub-project**: profile HMM
construction (deferred until a `locus_name` family has enough curated
diversity to support one), *de novo* gene structure prediction via
BRAKER/AUGUSTUS/Helixer (sub-project 3), ML/LLM-based sequence
refinement (sub-project 4), and RNA-seq integration.

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
which already has 12 examples across Pezizomycotina/Mucoromycota),
building real profile HMMs for those specific families becomes a
natural, separate follow-up — not a v1 blocker.

## Pipeline stages

### 1. Taxonomy routing (optional)

If a taxid is supplied, resolve its lineage via taxonkit (reusing
`MATPredict.db.taxonomy.resolve_lineage`) and match it against
`db/<Phylum>/order.yml`'s locus_name entries for that phylum, narrowing
the search to only the curated proteins from that phylum's
architecture(s). If no taxid is supplied, or its lineage doesn't
resolve to a phylum with curated coverage, fall back to searching the
full curated protein set across all phyla and let whichever family
scores best win -- slower, but keeps the tool usable on
unclassified/novel organisms.

### 2. Search

**Fast path (existing gene predictions available)**: search the
genome's predicted proteins directly against the curated reference
protein set (diamond or blastp). This is the common case when a
draft annotation already exists and is the cheaper, more accurate path
when available.

**Fallback path (genome FASTA only, no predictions)**: use
spliced protein-to-genome alignment -- exonerate
(`--model protein2genome`) or minimap2 (`-x splice`) -- to align each
curated reference protein directly against the genomic sequence.
Naive six-frame translation was considered and rejected: curated
records already show many intron-containing MAT genes (e.g.
Cryptococcus SXI1, Coprinopsis HD1), which six-frame translation would
fragment or miss entirely.

**Dual-anchor targets**: both `core_MAT` and `flanking_conserved` genes
from the matched reference record(s) are searched independently in the
same pass -- neither is a prerequisite for searching the other. This
matters because core MAT genes are frequently more divergent
(harder to detect at standard homology thresholds) than their flanking
genes, which are often broadly conserved across a wider taxonomic
range.

### 3. Boundary calling and confidence tiering

- If both `core_MAT` and `flanking_conserved` hits are found within a
  window sized off the matched reference record's own `locus.core`
  span (or a generous multiple of it, when spans vary a lot within a
  family) on a consistent contig/strand: report the spanning region as
  the locus boundary, confidence **high**.
- If `flanking_conserved` hits are found and cluster together, but no
  `core_MAT` hit clears the standard threshold: do **not** conclude the
  core gene is absent. Run a second, more sensitive search restricted
  to the genomic window between/around the flanking hits specifically
  for the core_MAT gene (relaxed e-value threshold and/or exonerate's
  more permissive alignment modes). If this second pass finds a
  plausible hit, report it with the boundary from the flanking hits,
  confidence **medium**, and flag the core gene's identity as
  needing manual confirmation. If the second pass still finds nothing,
  report the flanking-defined region as a candidate locus with no
  confirmed core gene, confidence **medium**, core gene marked absent
  rather than guessed.
- If only a `core_MAT` hit is found with no supporting flanking hits
  nearby: report it, confidence **low**, no defined boundary (matches
  the original design's core-anchored fallback for cases with sparse
  or absent flanking-gene homology in very divergent lineages).

### 4. Output

- **GFF3**: predicted locus region and gene features, using the same
  1-based fully-closed coordinate convention and `role`/`present`
  semantics as sub-project 1's schema, so predictions and curated
  records are directly comparable.
- **Validation-style report**: per-gene identity/coverage against the
  matched reference protein(s), which specific curated record(s) the
  match was scored against, confidence tier, and which expected genes
  (per the matched `locus_name`'s `order.yml` entry) were not found --
  deliberately mirrors `metadata.yaml`'s `validation` block shape from
  sub-project 1, so a detection result can be read the same way a
  curated record's validation output is read.

## Sn/Sp benchmark suite (sub-project 6)

**Method**: leave-one-out against sub-project 1's own curated records.
For each accepted record with a real, fetchable source sequence
(`coordinate_provenance` other than `not_available`), temporarily
exclude that record from the reference protein set used for search,
then run `matpredict detect` against that record's own source
region (fetched fresh, not read from the curated GFF3/GBK, to avoid
circularity), and score whether the pipeline recovers:
- the correct gene identities (matches `genes[]` with `present: true`)
- the correct locus boundary (within some tolerance, compared to
  `locus.core`)
- the correct confidence tier (a record with a well-conserved locus
  should score high; one known to have an unusually divergent core
  gene should reasonably score medium)

Records marked `excluded_from_coordinate_benchmark: true` (the Z.
rouxii translocation-background record, the Taphrina tblastn-only
record) are excluded from boundary-accuracy scoring specifically, but
remain usable for gene-identity-only scoring.

**Reporting**: Sn/Sp computed per `locus_name` family (not pooled
across the whole DB), since detection difficulty genuinely differs by
architecture -- a family with 12 examples (`MAT`) is a much fairer test
of the search method than a family with 1 (`MATyl`), and pooling would
hide that.

## Testing / acceptance criteria for this sub-project

- Unit tests for taxonomy routing, boundary-calling logic (all three
  confidence tiers, including the flanking-only-with-relaxed-core-search
  case), and GFF3/report output formatting -- using mocked search-tool
  output, not live diamond/exonerate runs, matching sub-project 1's
  established testing pattern (never call live external tools from unit
  tests).
- The leave-one-out benchmark itself, run against the real curated DB,
  as an integration-level check (not a fast unit test) -- reports Sn/Sp
  per family, not a pass/fail gate, since some families are expected to
  score poorly with only 1-2 examples until more curation happens.
- Acceptance target: the pipeline runs end-to-end (taxonomy routing →
  search → boundary calling → GFF3 + report output) against at least
  one real held-out record from each of the three phyla, and the
  leave-one-out benchmark produces a real, non-fabricated Sn/Sp report
  covering every family with 2+ examples.

## Out of scope for this sub-project

- Profile HMM construction for any `locus_name` family (a natural
  follow-up once specific families accumulate enough diversity)
- *De novo* gene structure prediction (BRAKER/AUGUSTUS/Helixer) --
  sub-project 3
- ML/LLM-based boundary/motif refinement -- sub-project 4
- RNA-seq integration
- A public web resource for browsing detection results -- sub-project 5
