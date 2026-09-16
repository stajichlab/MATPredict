# MATPredict Sub-project 1: Curated MAT Locus Reference Database

Date: 2026-09-16 (revised after Opus review)
Status: draft, pending user review

## Context

MATPredict aims to build a general-purpose MAT (mating-type) locus
annotator/classifier for fungi. MAT locus structure varies enough
across phyla that detection, gene prediction, and classification all
need clade-specific models. All downstream work (HMMs, synteny models,
gene prediction refinement, ML/LLM classification, a public web
resource, and Sn/Sp benchmarking) depends on a curated, validated
ground-truth set of MAT loci. This is too large a project for one
spec, so it is decomposed into sub-projects:

1. **Reference database** (this spec)
2. Detection/annotation pipeline (homology + synteny, taxonomy-aware)
3. Gene prediction refinement (clade-specific trained gene models)
4. ML/LLM classifier layer (boundary refinement, motif confirmation, idiomorph classification)
5. Web resource / publication
6. Validation/benchmark suite (Sn/Sp), designed alongside 1 and 2

This spec covers only sub-project 1.

Existing repo state this builds on: `db/{Ascomycota,Basidiomycota,Mucoromycota}/`
currently holds rough `annotated.yml` (per-organism gene/accession lists)
and `order.yml` (per-clade expected gene complement per mating type)
stubs. `testset/Zygo/` has a working prototype for one clade: homology
search (diamond) against known MAT proteins, presence/absence calling
(`abspres.csv`), followed by alignment/tree building to validate
orthology. `testset/get_MAT.py` is an unfiltered NCBI keyword scraper
(searches "mating-type", downloads without curation) — not part of
this design; it predates the curated-DB approach and is superseded by
the workflow below.

**Migration note**: `annotated.yml` becomes the seed source for the
first round of literature-mining candidates (its accessions/organisms
get re-verified through the standard workflow, not copied in
directly). `order.yml`'s per-clade expected gene complement becomes the
authoritative source for the `mating_type.idiomorph` and `genes[].role`
controlled vocabulary per phylum — referenced from, not duplicated
into, `_schema/metadata.schema.yaml` (the schema loads clade-specific
enums from `db/<Phylum>/order.yml` at validation time). Both legacy
files are left in place until every organism they reference has been
re-curated into the new layout or explicitly rejected, then removed.

## Scope

Populate a curated MAT locus reference database across all three
phyla in parallel (Ascomycota, Basidiomycota, Mucoromycota), restricted
initially to **tier-1 evidence**: loci from peer-reviewed publications
with experimental validation, not computational inference alone.
Tier-2 (homology-inferred) candidates are tracked in the schema but not
admitted to the accepted database in this phase. Evidence tiering is
per-claim, not per-record (see schema below), since a record's locus
existence, boundaries, and idiomorph assignment often carry different
evidence strength even within one tier-1 publication.

**Done criteria for this sub-project**: schema, curation tooling, and
validation pipeline built and proven on 5-10 tier-1 records per phylum
(~15-30 records total), each passing all validation checks below,
deliberately including at least one coordinate-less (genetic/RFLP-only)
record and one pre-genomic INSDC-nucleotide-only record per phylum
where literature allows, so the schema's edge-case paths are exercised
before scaling up.

## Tool architecture (applies to this and all later sub-projects)

MATPredict converges on a single installable Python package with one
CLI and subcommands per stage (`matpredict curate-db`, `matpredict
detect`, `matpredict predict`, `matpredict classify`), matching the
existing `src/MATPredict/__main__.py` stub. Shared infrastructure
(taxonomy resolution via taxonkit, schema validation, GFF3/GBK/YAML
I/O) lives in one place under `src/MATPredict/`. A menu-driven
interactive mode is itself a subcommand, not a separate tool. Rust is
not a parallel track — it enters later only as an optional
subprocess-called accelerator for a specific hot path proven too slow
in Python (e.g. a synteny scan), once training/parameters are
codified.

This sub-project's curation and validation code is built as real
modules — `src/MATPredict/db/schema.py`, `src/MATPredict/db/curate.py`,
`src/MATPredict/db/validate.py`, `src/MATPredict/db/build_duckdb.py` —
wired to a `matpredict curate-db` subcommand, not as standalone
one-off scripts.

## Directory layout

```
db/
  <Phylum>/
    order.yml                    # per-clade idiomorph + gene-role controlled vocabulary (existing file, retained)
    <Order_or_Family>/
      <taxid>_<strain_slug>_<idiomorph>/
        locus.gff3                # accepted records only; core MAT locus + flanking genes
        locus.gbk                 # accepted records only; full annotated GenBank record for the locus region
        proteins.faa               # accepted records only; extracted protein FASTA for each gene in genes[] (sub-project 2's actual input)
        metadata.yaml              # provenance, evidence, taxonomy, validation, versioning
  candidates/
    <Phylum>/
      <proposed_record_id>/
        metadata.yaml              # status: needs_review | rejected; no gff3/gbk/proteins.faa until accepted
  _schema/
    metadata.schema.yaml           # schema metadata.yaml must satisfy
    duckdb_schema.sql              # DDL for the derived query database
  _release.yml                     # release_tag -> [record_id@record_version, ...] manifest
```

Record key: `taxid + strain/isolate + idiomorph`. A heterothallic
species with two idiomorphs gets two records; a homothallic species
with both idiomorphs present gets one record spanning both, per how
the source publication describes the locus. `record_id` is immutable
once assigned (even if the strain name is later corrected) — corrections
go in `organism.strain.name` with the original slug preserved in the
directory path and record_id.

`needs_review`/`rejected` candidates live under `db/candidates/`, not
in the accepted phylum tree, so the accepted `db/<Phylum>/` tree only
ever contains records that passed the human review gate.

## `metadata.yaml` schema

```yaml
record_id: <taxid>_<strain_slug>_<idiomorph>     # matches directory name; immutable once assigned
record_version: 1              # increments on any edit to core/extended_flank/genes/evidence

taxonomy:
  taxid: 4837
  lineage: "k__Fungi;p__Mucoromycota;...;s__Phycomyces_blakesleeanus"  # taxonkit-resolved, cached
  lineage_resolved_date: 2026-09-16

organism:
  species: "Phycomyces blakesleeanus"
  strain:
    name: "NRRL 1555"
    known: true                        # false when the paper doesn't name a strain
    culture_collection_ids: ["NRRL 1555", "CBS 253.65"]  # cross-referenced synonyms, if any
    differs_from_sequenced: false      # true if the published strain != the strain the cited assembly represents

mating_type:
  idiomorph: "Plus"          # controlled vocabulary loaded from db/<Phylum>/order.yml
  system: "heterothallic"    # heterothallic | homothallic | pseudohomothallic

locus:
  coordinate_provenance: "published_explicit"  # published_explicit | curator_derived | not_available
  # not_available records (e.g. RFLP/genetic-mapping-only evidence) omit `core` entirely and are
  # flagged excluded_from_coordinate_benchmark: true, so sub-project 6 knows not to score them on
  # boundary Sn/Sp, while they remain usable for gene-content/citation purposes.
  excluded_from_coordinate_benchmark: false
  core:                       # boundary as defined by the source publication; never overwritten in place
    completeness: "complete"   # complete | partial | fragmented
    reference_orientation: "tptA->rnhA"  # defined by cited flanking-gene order, not a raw +/- strand call,
                                          # since idiomorphs at one locus are often on opposite strands
    definition_note: "core boundary per Idnurm et al. 2008, between tptA and rnhA"
    segments:                  # one entry normally; >1 when the locus is split across contigs
      - sequence_source:
          type: "assembly"       # assembly | insdc_nucleotide | none
          accession: "GCA_000315115.1"   # fully versioned
          seq_region: "scaffold_3"
        start: 120345           # 1-based, fully-closed, matching GFF3 convention (documented once, enforced by validator)
        end: 128900
        contig_edge_distance: null   # bp from segment end to contig/scaffold end; null if not near an edge
        sequence_checksum: "md5:1a2b3c..."  # md5 of the extracted segment sequence, computed at curation time
  extended_flank: null        # populated later when synteny work walks boundaries outward; additive, versioned

genes:
  - gene_index: 0               # stable ordinal within this record; PK component (handles paralogs/duplicates)
    name: sexP
    protein_accession: "A0A078N0N5_9FUNG"
    role: "core_MAT"            # core_MAT | flanking_conserved | flanking_variable, loaded from order.yml per clade
    present: true                # false = curated absence call (needed as true negatives for sub-project 6)
    locus_tag: null
    start: 121002
    end: 122400
    strand: "+"
    order_in_locus: 1
  - gene_index: 1
    name: tptA
    protein_accession: "B0F2G7_PHYBL"
    role: "flanking_conserved"
    present: true
    locus_tag: null
    start: 120345
    end: 121000
    strand: "+"
    order_in_locus: 0

evidence:
  # tiered per claim, not per record: a tier-1 paper can still leave one of these at tier 2/uncertain
  locus_existence:
    tier: 1
    citations:
      - pmid: "18248337"
        doi: "10.1128/EC.00281-08"
    experimental_method: "targeted sequencing + genetic crosses"
  boundaries:
    tier: 1
    citations:
      - pmid: "18248337"
  idiomorph_assignment:
    tier: 1
    citations:
      - pmid: "18248337"

validation:
  accession_resolved: true
  accession_resolved_date: 2026-09-16
  accession_resolved_version: "GCA_000315115.1"   # exact version that passed, not just a boolean+date
  sequence_match:
    status: "pass"              # pass | warn | fail (tri-state, not boolean)
    percent_identity: 99.8
    coverage: 100.0
    notes: ""
  taxonomy_current: true
  status: "accepted"            # accepted | needs_review | rejected
  rejection_reason: null        # required non-null when status == "rejected"

curation:
  proposed_by: "literature-mining-agent"
  reviewed_by: "jason.stajich@ucr.edu"
  reviewed_date: 2026-09-16

model_provenance: null         # reserved: populated by sub-project 2/4 when a model is trained citing this record
```

Fields not yet known at curation time (e.g. an accession the paper
doesn't specify) are left blank, never inferred or guessed.

**Coordinate convention**: all positions are 1-based, fully-closed
(GFF3 convention), asserted by the schema validator on every load.

## Release manifest (replaces per-file `db_release` stamping)

`db/_release.yml`:

```yaml
releases:
  "2026.09.0":
    cut_date: 2026-09-16
    records:
      - "4837_NRRL-1555_Plus@1"
      - "4837_NRRL-1555_Minus@1"
      # ...
```

Cutting a release appends one entry here rather than rewriting
`record_version`/`db_release` into every `metadata.yaml`. A model
(sub-project 2/4) cites `model_provenance.db_release: "2026.09.0"` and
that manifest entry pins the exact `record_id@record_version` set it
was trained on, even as later releases add or revise records.

## DuckDB schema (query cache, not source of truth)

Generated by a build script that walks all `metadata.yaml` files
(both `db/<Phylum>/` accepted records and `db/candidates/`).
GFF3/GBK/YAML/`_release.yml` on disk remain authoritative; DuckDB is
rebuilt from them, never edited directly. The build script — not
DuckDB itself — enforces referential integrity across tables.
`phylum`/`order_or_family` are derived from `taxonomy.lineage` and the
file path at build time (never hand-edited, so they can't drift from
the source of truth independently).

```sql
CREATE TABLE organism (
    taxid INTEGER PRIMARY KEY,
    species TEXT NOT NULL,
    lineage TEXT NOT NULL,
    lineage_resolved_date DATE
);

CREATE TABLE locus_record (
    record_id TEXT PRIMARY KEY,
    taxid INTEGER REFERENCES organism(taxid),
    strain_name TEXT,
    strain_known BOOLEAN,
    idiomorph TEXT NOT NULL,
    mating_system TEXT,                -- heterothallic | homothallic | pseudohomothallic
    phylum TEXT NOT NULL,              -- derived from lineage at build time
    order_or_family TEXT,              -- derived from directory path at build time
    record_version INTEGER NOT NULL,
    coordinate_provenance TEXT NOT NULL,   -- published_explicit | curator_derived | not_available
    excluded_from_coordinate_benchmark BOOLEAN NOT NULL DEFAULT FALSE,
    completeness TEXT,                 -- complete | partial | fragmented
    reference_orientation TEXT,
    definition_note TEXT,
    evidence_tier_locus_existence INTEGER,
    evidence_tier_boundaries INTEGER,
    evidence_tier_idiomorph INTEGER,
    experimental_method TEXT,
    validation_status TEXT NOT NULL,   -- accepted | needs_review | rejected
    rejection_reason TEXT,
    accession_resolved BOOLEAN,
    accession_resolved_version TEXT,
    sequence_match_status TEXT,        -- pass | warn | fail
    sequence_match_pct_identity DOUBLE,
    sequence_match_coverage DOUBLE,
    taxonomy_current BOOLEAN,
    proposed_by TEXT,
    reviewed_by TEXT,
    reviewed_date DATE,
    gff3_path TEXT,                    -- null for candidates (not yet generated)
    gbk_path TEXT,
    proteins_fasta_path TEXT,
    metadata_path TEXT NOT NULL
);
CREATE INDEX idx_locus_record_taxid ON locus_record(taxid);
CREATE INDEX idx_locus_record_phylum ON locus_record(phylum);
CREATE INDEX idx_locus_record_idiomorph ON locus_record(idiomorph);
CREATE INDEX idx_locus_record_status ON locus_record(validation_status);

CREATE TABLE locus_segment (
    record_id TEXT REFERENCES locus_record(record_id),
    segment_index INTEGER NOT NULL,
    sequence_source_type TEXT NOT NULL,   -- assembly | insdc_nucleotide | none
    accession TEXT,                        -- fully versioned
    seq_region TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    contig_edge_distance INTEGER,
    sequence_checksum TEXT,
    PRIMARY KEY (record_id, segment_index)
);

CREATE TABLE locus_gene (
    record_id TEXT REFERENCES locus_record(record_id),
    gene_index INTEGER NOT NULL,       -- handles paralogs/duplicate gene names within one record
    name TEXT NOT NULL,
    protein_accession TEXT,
    role TEXT NOT NULL,                -- core_MAT | flanking_conserved | flanking_variable
    present BOOLEAN NOT NULL,          -- false = curated absence call (true negative for sub-project 6)
    locus_tag TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    strand TEXT,
    order_in_locus INTEGER,
    PRIMARY KEY (record_id, gene_index)
);

CREATE TABLE citation (
    record_id TEXT REFERENCES locus_record(record_id),
    claim TEXT NOT NULL,               -- locus_existence | boundaries | idiomorph_assignment
    pmid TEXT,
    doi TEXT,
    PRIMARY KEY (record_id, claim, pmid)
);

CREATE TABLE extended_flank (
    record_id TEXT REFERENCES locus_record(record_id),
    walked_out_date DATE,
    seq_region TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    method TEXT,                       -- how the extension was derived (synteny tool, manual, etc.)
    is_current BOOLEAN NOT NULL DEFAULT TRUE,  -- superseded extensions kept with is_current=FALSE, not deleted
    notes TEXT,
    PRIMARY KEY (record_id, walked_out_date)
);

CREATE TABLE release_manifest (
    release_tag TEXT NOT NULL,
    record_id TEXT NOT NULL,
    record_version INTEGER NOT NULL,
    cut_date DATE,
    PRIMARY KEY (release_tag, record_id)
);
```

## Curation workflow

Enforces verifiability at every step; no record reaches `accepted`
without passing automated validation and explicit human sign-off.

1. **Literature search** (per phylum/order, via PubMed, seeded in the
   first pass by re-verifying `annotated.yml` entries): agent proposes
   candidates under `db/candidates/<Phylum>/`. Every candidate must
   carry a resolvable PMID/DOI and the specific text/table the claim
   was drawn from. No candidate is proposed without one. Candidates
   with genetic/RFLP-only evidence and no sequence coordinates are
   proposed with `coordinate_provenance: not_available` rather than
   discarded.
2. **Candidate extraction** → draft `metadata.yaml` with
   `validation.status="needs_review"` and per-claim evidence tiers set
   from what the paper actually supports. Fields not explicitly stated
   in the source are left blank, never inferred.
3. **Automated validation** (NCBI E-utilities + taxonkit), skipped for
   fields that don't apply to a given `coordinate_provenance`:
   - `accession_resolved` (+ `accession_resolved_version`): does the
     cited accession still resolve, at which version, and is it not
     suppressed/replaced?
   - `sequence_match`: percent identity + coverage between the
     assembly region's translation and the cited protein sequence,
     scored pass/warn/fail against thresholds (not a strict boolean,
     since old GenBank translations and re-called ORFs commonly
     diverge slightly).
   - `taxonomy_current`: does taxonkit resolve the species to a
     current, non-merged taxid?
   - Any `fail` or unresolved check keeps `validation.status =
     "needs_review"`; no auto-accept. `warn` is surfaced to the
     reviewer but doesn't block review.
4. **Human review gate** (user): review `needs_review` records in
   `db/candidates/` with citations + validation results side by side.
   Accept — move the record directory from `db/candidates/<Phylum>/`
   to `db/<Phylum>/<Order_or_Family>/`, set `status="accepted"` — or
   reject — set `status="rejected"` with a required
   `rejection_reason`, leave it under `db/candidates/` for traceability
   so the mining agent doesn't re-propose it blindly.
5. **GFF3/GBK/protein FASTA generation**: only for accepted records,
   built from the validated coordinates/accessions. `proteins.faa` is
   the actual per-record input sub-project 2 consumes.
6. **Release cut** (periodic, manual trigger): append the current set
   of `accepted` `record_id@record_version` pairs to `db/_release.yml`
   under a new release tag.

## Testing / acceptance criteria

- **Schema validation**: every `metadata.yaml` (accepted and
  candidate) validates against `metadata.schema.yaml`, with clade
  controlled vocabulary loaded from the matching `db/<Phylum>/order.yml`.
- **Cross-file consistency**: for `coordinate_provenance:
  published_explicit`/`curator_derived` records, `locus.gff3`
  coordinates match `metadata.yaml locus.core.segments`; gene
  coordinates/names in `locus.gff3`/`locus.gbk`/`proteins.faa` match
  `genes[]` entries with `present: true`.
- **Edge-case coverage**: at least one `not_available`-coordinate
  record, one `insdc_nucleotide`-sourced record, and one record with
  `present: false` genes exist in the seed set and pass their
  applicable checks (not_available records skip the coordinate
  cross-file check by design).
- **DuckDB build**: build script loads all records (accepted +
  candidates) without error; row counts match file counts; foreign
  keys the build script enforces are validated with no orphans.
- **Validation script unit tests**: mock NCBI responses to test
  accession-resolved / sequence-match (pass/warn/fail thresholds) /
  taxonomy-current logic independently of live network calls.
- **Acceptance target**: 5-10 tier-1 records per phylum, including the
  edge cases above, all passing the checks, committed under
  `db/<Phylum>/`.

## Out of scope for this sub-project

- Synteny-based boundary walking (`extended_flank` population) —
  schema supports it, but population is deferred to sub-project 2.
- Tier-2 (homology-only) candidate acceptance.
- Model training / `model_provenance` population — sub-project 2/4.
- Web resource / public browsing UI — sub-project 5.
- `testset/get_MAT.py` unfiltered scraper — superseded by the
  literature-mining workflow above; not modified in this sub-project.
- Full removal of `annotated.yml`/`order.yml` — `order.yml` is
  actively referenced (see above); both are removed only once fully
  migrated, as a follow-up cleanup, not part of this sub-project's
  acceptance criteria.
