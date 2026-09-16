# MATPredict Sub-project 1: Curated MAT Locus Reference Database

Date: 2026-09-16 (revised after Opus review, then Fable review)
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
currently holds a rough `annotated.yml` (per-organism gene/accession
lists) and a partial `order.yml` (per-clade expected gene complement
per mating type) for Mucoromycota only; Ascomycota and Basidiomycota
have no `order.yml` yet. `testset/Zygo/` has a working prototype for
one clade: homology search (diamond) against known MAT proteins,
presence/absence calling (`abspres.csv`), followed by alignment/tree
building to validate orthology. `testset/get_MAT.py` is an unfiltered
NCBI keyword scraper — not part of this design; superseded by the
workflow below. `src/MATPredict/__main__.py` is a copied CLI template
that does not currently run (`main()` dispatches to an undefined
`_menu_map_reads`); there is no `pyproject.toml`/`pixi.toml` yet.

**Migration note**: `annotated.yml` becomes the seed source for the
first round of literature-mining candidates (re-verified through the
standard workflow, not copied in directly). `order.yml` becomes the
authoritative source for each phylum's locus/idiomorph controlled
vocabulary (see "Locus and idiomorph model" below) and must be written
(or rewritten, for Mucoromycota) against a defined schema as part of
this sub-project's scope — see Scope. Both legacy files are removed
only once every organism they reference has been re-curated or
explicitly rejected.

## Scope

Populate a curated MAT locus reference database across all three
phyla in parallel (Ascomycota, Basidiomycota, Mucoromycota), restricted
initially to **tier-1 evidence**: loci from peer-reviewed publications
with experimental validation, not computational inference alone.
Tier-2 (homology-inferred) candidates are tracked in the schema but not
admitted to the accepted database in this phase. Evidence tiering is
per-claim, not per-record (see schema below).

This sub-project's scope explicitly includes, as prerequisite
implementation tasks (not deferred, not assumed to already exist):

- Fixing `src/MATPredict/__main__.py` into a working `matpredict`
  entry point with real subcommand dispatch (starting with
  `matpredict curate-db`), plus a `pyproject.toml`/`pixi.toml` so the
  package installs and runs.
- A small config layer: where the `db/` root lives (repo-relative, not
  package data), NCBI E-utilities email/API key, and an on-disk HTTP
  response cache for E-utilities calls (curation reruns should not
  re-hit NCBI for records already resolved).
- Writing `db/<Phylum>/order.yml` for all three phyla against a new
  `_schema/order.schema.yaml`, including rewriting Mucoromycota's
  existing file to the new schema (see "Locus and idiomorph model").

**Done criteria for this sub-project**: schema, curation tooling, and
validation pipeline built and proven on 5-10 tier-1 records per phylum
(~15-30 records total), each passing all validation checks below,
deliberately including per phylum, where literature allows: one
coordinate-less (genetic/RFLP-only) record, one pre-genomic
INSDC-nucleotide-only record, one record with `present: false` genes,
and for Basidiomycota specifically, one tetrapolar species with both
an `HD` and a `PR` locus record for the same strain — so the schema's
edge-case paths, including the multiallelic locus model, are exercised
before scaling up.

## Tool architecture (applies to this and all later sub-projects)

MATPredict converges on a single installable Python package with one
CLI and subcommands per stage (`matpredict curate-db`, `matpredict
detect`, `matpredict predict`, `matpredict classify`). Shared
infrastructure (taxonomy resolution via taxonkit, schema validation,
GFF3/GBK/YAML I/O, the config layer and HTTP cache above) lives in one
place under `src/MATPredict/`. A menu-driven interactive mode is
itself a subcommand, not a separate tool. Rust is not a parallel
track — it enters later only as an optional subprocess-called
accelerator for a specific hot path proven too slow in Python, once
training/parameters are codified.

This sub-project's curation and validation code is built as real
modules — `src/MATPredict/db/schema.py`, `src/MATPredict/db/curate.py`,
`src/MATPredict/db/validate.py`, `src/MATPredict/db/build_duckdb.py`,
`src/MATPredict/config.py` — wired to the `matpredict curate-db`
subcommand.

## Locus and idiomorph model

Not every clade has one locus with a two-way idiomorph. This schema
must represent three cases:

- **Bipolar heterothallic** (most Ascomycota, e.g. MAT1-1/MAT1-2; many
  Mucoromycota, Plus/Minus): one locus, one of a small closed set of
  idiomorph values.
- **Homothallic**: one locus, both idiomorphs present together in one
  genome region.
- **Tetrapolar Basidiomycota**: two physically unlinked loci per
  strain (conventionally named `HD` and `PR`), each independently
  **multiallelic** (A1, A2, ... An — an open-ended set, not a fixed
  enum).

To cover all three: `record_id` includes a `locus_name` component
(default `"MAT"` for the single-locus bipolar/homothallic case; `"HD"`
or `"PR"` for tetrapolar Basidiomycota), and `mating_type.idiomorphs`
is always a list (one element for heterothallic, two for a combined
homothallic record). `db/<Phylum>/order.yml` declares, per
`locus_name`, whether idiomorph values are a closed enum
(`vocabulary_type: enum`, with the allowed list) or an open,
pattern-validated allele set (`vocabulary_type: pattern`, e.g.
`^A[0-9]+$`) — the schema validator enforces whichever `order.yml`
declares for that phylum/locus_name, never a single hardcoded
Plus/Minus/MAT1-1/MAT1-2 enum baked into the metadata schema itself.

`db/<Phylum>/order.yml` (new schema, validated by
`_schema/order.schema.yaml`):

```yaml
phylum: Mucoromycota
loci:
  - locus_name: "MAT"
    vocabulary_type: "enum"
    idiomorph_values: ["Plus", "Minus"]
    genes:
      - name: tptA
        role: flanking_conserved
      - name: sexP
        role: core_MAT
        present_in_idiomorphs: ["Plus"]
      - name: sexM
        role: core_MAT
        present_in_idiomorphs: ["Minus"]
      - name: rnhA
        role: flanking_conserved
```

For a tetrapolar Basidiomycota phylum, two `loci` entries (`HD`, `PR`)
each with `vocabulary_type: pattern` and a `idiomorph_pattern` instead
of `idiomorph_values`.

## Directory layout

```
db/
  <Phylum>/
    order.yml                    # per-clade locus/idiomorph/gene-role controlled vocabulary
    <Order_or_Family>/
      <taxid>_<strain_slug>_<locus_name>_<idiomorph_key>/
        locus.gff3                # accepted records only
        locus.gbk                 # accepted records only
        proteins.faa               # accepted records only; header convention below
        metadata.yaml
  candidates/
    <Phylum>/
      <proposed_record_id>/
        metadata.yaml              # status: needs_review | rejected; no gff3/gbk/proteins.faa until accepted
  _schema/
    metadata.schema.yaml
    order.schema.yaml
    duckdb_schema.sql
  _release.yml
```

`idiomorph_key` is the path-safe form of `mating_type.idiomorphs`:
the single value for a heterothallic record (e.g. `Plus`, `A1`), or
`combined` for a homothallic record spanning both idiomorphs (the full
list is still recorded in `mating_type.idiomorphs` in the metadata).

`record_id` is immutable once assigned. Editing an already-accepted
record's `locus.core`, `genes[]`, or any `evidence.*.tier` is a
**material change**: it bumps `record_version` in place (record stays
under `db/<Phylum>/`, full history via git) and requires the same
human review gate as a new candidate before the new version's
`validation.status` can read `accepted` again — it is not routed back
through `db/candidates/`. Metadata-only corrections (e.g. fixing a
typo in `organism.strain.name`) bump `record_version` without
re-review. The validator enforces the invariant that
`validation.status == "accepted"` if and only if the record's
directory is under `db/<Phylum>/...` (never under `db/candidates/`),
so status and location can never disagree.

Strain slugging rule: lowercase, spaces and `/` become `-`, other
punctuation stripped; `known: false` strains slug to `unknown-<n>`
where `<n>` is a per-species counter, so multiple unnamed-strain
records for one species don't collide.

## `metadata.yaml` schema

```yaml
record_id: <taxid>_<strain_slug>_<locus_name>_<idiomorph_key>   # matches directory name; immutable once assigned
record_version: 1              # increments on any edit; material changes require re-review to reach "accepted" again

taxonomy:
  taxid: 4837
  lineage: "k__Fungi;p__Mucoromycota;...;s__Phycomyces_blakesleeanus"  # taxonkit-resolved, cached
  lineage_resolved_date: 2026-09-16

organism:
  species: "Phycomyces blakesleeanus"
  strain:
    name: "NRRL 1555"
    known: true                        # false when the paper doesn't name a strain
    culture_collection_ids: ["NRRL 1555", "CBS 253.65"]
    differs_from_sequenced: false      # true if published strain != the strain the cited assembly represents

mating_type:
  locus_name: "MAT"                    # "MAT" (single-locus) | "HD" | "PR" (tetrapolar Basidiomycota)
  idiomorphs: ["Plus"]                 # list: 1 element (heterothallic/tetrapolar), 2 (homothallic combined)
  system: "heterothallic"              # heterothallic | homothallic | pseudohomothallic
  # idiomorph value format (enum member vs. pattern-matched allele id) is validated against
  # the matching locus_name entry in db/<Phylum>/order.yml, not a hardcoded schema enum.

locus:
  coordinate_provenance: "published_explicit"  # published_explicit | curator_derived | not_available
  excluded_from_coordinate_benchmark: false     # true whenever coordinate_provenance != published_explicit-with-passing-checks
  core:
    completeness: "complete"   # complete | partial | fragmented
    reference_orientation: "tptA->rnhA"  # defined by cited flanking-gene order, not a raw +/- strand call
    definition_note: "core boundary per Idnurm et al. 2008, between tptA and rnhA"
    segments:                  # one entry normally; >1 when the locus is split across contigs
      - segment_index: 0
        sequence_source:
          type: "assembly"       # assembly | insdc_nucleotide | none
          accession: "GCA_000315115.1"   # fully versioned; ncbi_protein:/uniprotkb: prefix convention for gene accessions below
          seq_region: "scaffold_3"
        start: 120345           # 1-based, fully-closed (GFF3 convention), enforced by validator
        end: 128900
        contig_edge_distance: null
        sequence_checksum: "md5:1a2b3c..."
  extended_flank: []           # list, not a nullable scalar; empty until synteny work (sub-project 2) appends an entry;
                                # superseded entries are kept with is_current: false, never deleted

genes:
  - gene_index: 0
    name: sexP
    protein_accession: "ncbi_protein:AAB12345.1"   # namespaced: ncbi_protein: | uniprotkb:, resolver chosen by prefix
    role: "core_MAT"            # validated against db/<Phylum>/order.yml for this locus_name
    present: true                # false = curated absence call (true negative for sub-project 6)
    locus_tag: null
    segment_index: 0             # which locus.core.segments entry this gene falls on
    start: 121002
    end: 122400
    strand: "+"
    order_in_locus: 1
  - gene_index: 1
    name: tptA
    protein_accession: "ncbi_protein:CAB67890.1"
    role: "flanking_conserved"
    present: true
    locus_tag: null
    segment_index: 0
    start: 120345
    end: 121000
    strand: "+"
    order_in_locus: 0

evidence:
  # tiered and cited per claim, not per record
  locus_existence:
    tier: 1
    experimental_method: "targeted sequencing + genetic crosses"
    citations: [{pmid: "18248337", doi: "10.1128/EC.00281-08"}]
  boundaries:
    tier: 1
    experimental_method: "targeted sequencing"
    citations: [{pmid: "18248337", doi: null}]
  idiomorph_assignment:
    tier: 1
    experimental_method: "genetic crosses"
    citations: [{pmid: "18248337", doi: null}]

validation:
  accession_resolved: true
  accession_resolved_date: 2026-09-16
  accession_resolved_version: "GCA_000315115.1"
  sequence_match:
    status: "pass"              # pass | warn | fail; record-level status = worst of all per-gene checks below
    per_gene:
      - gene_index: 0
        percent_identity: 99.8
        coverage: 100.0
        status: "pass"
    notes: ""
  taxonomy_current: true
  status: "accepted"            # accepted | needs_review | rejected
  rejection_reason: null        # required non-null when status == "rejected"

curation:
  proposed_by: "literature-mining-agent"
  proposal_dedupe_key: "18248337|4837|MAT"   # pmid|taxid|locus_name, used to suppress re-proposing a rejected candidate
  reviewed_by: "jason.stajich@ucr.edu"
  reviewed_date: 2026-09-16

model_provenance: null         # reserved: populated by sub-project 2/4 when a model is trained citing this record
```

Fields not yet known at curation time are left blank, never inferred
or guessed. **Coordinate convention**: all positions are 1-based,
fully-closed (GFF3 convention), asserted by the validator on every
load.

**`proteins.faa` header convention**: one entry per `genes[]` row with
`present: true`, header
`>{record_id}|gene_index={gene_index}|name={name}|role={role}` — this
is the literal input sub-project 2 parses.

## Release manifest

`db/_release.yml`:

```yaml
releases:
  "2026.09.0":
    cut_date: 2026-09-16
    git_tag: "db-release-2026.09.0"    # required; `git tag` is cut at the same commit, making the release's exact
                                        # file content recoverable, not just the record_id@version list
    records:
      - "4837_nrrl-1555_MAT_Plus@1"
      - "4837_nrrl-1555_MAT_Minus@1"
```

## DuckDB schema (query cache, not source of truth)

Generated by a build script that walks all `metadata.yaml` files
(`db/<Phylum>/` and `db/candidates/`). GFF3/GBK/YAML/`_release.yml` on
disk remain authoritative; DuckDB is rebuilt from them, never edited
directly. The build script enforces referential integrity (DuckDB does
not enforce FKs). `phylum`/`order_or_family` are derived from
`taxonomy.lineage` and the file path at build time, never hand-edited.

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
    locus_name TEXT NOT NULL,          -- MAT | HD | PR
    idiomorph_key TEXT NOT NULL,       -- path-safe form; full list is in locus_idiomorph
    mating_system TEXT,                -- heterothallic | homothallic | pseudohomothallic
    phylum TEXT NOT NULL,              -- derived from lineage at build time
    order_or_family TEXT,              -- derived from directory path at build time
    record_version INTEGER NOT NULL,
    coordinate_provenance TEXT NOT NULL,
    excluded_from_coordinate_benchmark BOOLEAN NOT NULL DEFAULT FALSE,
    completeness TEXT,
    reference_orientation TEXT,
    definition_note TEXT,
    validation_status TEXT NOT NULL,   -- accepted | needs_review | rejected
    rejection_reason TEXT,
    accession_resolved BOOLEAN,
    accession_resolved_version TEXT,
    sequence_match_status TEXT,        -- pass | warn | fail
    taxonomy_current BOOLEAN,
    proposed_by TEXT,
    proposal_dedupe_key TEXT,
    reviewed_by TEXT,
    reviewed_date DATE,
    gff3_path TEXT,                    -- null for candidates
    gbk_path TEXT,
    proteins_fasta_path TEXT,
    metadata_path TEXT NOT NULL
);
CREATE INDEX idx_locus_record_taxid ON locus_record(taxid);
CREATE INDEX idx_locus_record_phylum ON locus_record(phylum);
CREATE INDEX idx_locus_record_locus_name ON locus_record(locus_name);
CREATE INDEX idx_locus_record_status ON locus_record(validation_status);

CREATE TABLE locus_idiomorph (
    record_id TEXT REFERENCES locus_record(record_id),
    idiomorph_value TEXT NOT NULL,     -- one row per element of mating_type.idiomorphs
    PRIMARY KEY (record_id, idiomorph_value)
);

CREATE TABLE locus_segment (
    record_id TEXT REFERENCES locus_record(record_id),
    segment_index INTEGER NOT NULL,
    sequence_source_type TEXT NOT NULL,   -- assembly | insdc_nucleotide | none
    accession TEXT,
    seq_region TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    contig_edge_distance INTEGER,
    sequence_checksum TEXT,
    PRIMARY KEY (record_id, segment_index)
);

CREATE TABLE locus_gene (
    record_id TEXT REFERENCES locus_record(record_id),
    gene_index INTEGER NOT NULL,
    segment_index INTEGER,
    name TEXT NOT NULL,
    protein_accession TEXT,            -- namespaced: ncbi_protein:... | uniprotkb:...
    role TEXT NOT NULL,
    present BOOLEAN NOT NULL,
    locus_tag TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    strand TEXT,
    order_in_locus INTEGER,
    sequence_match_status TEXT,        -- per-gene pass | warn | fail
    percent_identity DOUBLE,
    coverage DOUBLE,
    PRIMARY KEY (record_id, gene_index)
);

CREATE TABLE evidence_claim (
    record_id TEXT REFERENCES locus_record(record_id),
    claim TEXT NOT NULL,               -- locus_existence | boundaries | idiomorph_assignment
    tier INTEGER NOT NULL,
    experimental_method TEXT,
    PRIMARY KEY (record_id, claim)
);

CREATE TABLE citation (
    citation_id TEXT PRIMARY KEY,      -- surrogate key: uuid or hash(record_id, claim, pmid, doi)
    record_id TEXT REFERENCES locus_record(record_id),
    claim TEXT NOT NULL,
    pmid TEXT,
    doi TEXT
);

CREATE TABLE extended_flank (
    flank_id TEXT PRIMARY KEY,         -- surrogate key; multiple same-day entries no longer collide
    record_id TEXT REFERENCES locus_record(record_id),
    walked_out_date DATE,
    seq_region TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    method TEXT,
    is_current BOOLEAN NOT NULL DEFAULT TRUE,
    notes TEXT
);

CREATE TABLE release_manifest (
    release_tag TEXT NOT NULL,
    git_tag TEXT NOT NULL,
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
   carry a resolvable PMID/DOI, the specific text/table the claim was
   drawn from, and a `proposal_dedupe_key` (`pmid|taxid|locus_name`)
   checked against existing `rejected` candidates before proposing, so
   a rejected candidate isn't blindly re-proposed. Candidates with
   genetic/RFLP-only evidence are proposed with
   `coordinate_provenance: not_available` rather than discarded.
2. **Candidate extraction** → draft `metadata.yaml` with
   `validation.status="needs_review"` and per-claim evidence tiers set
   from what the paper actually supports. Fields not explicitly stated
   are left blank, never inferred.
   - **2b. Coordinate derivation** (required before step 3 can run
     for anything other than `not_available`): if the paper cites a
     protein/gene accession but not explicit locus coordinates, the
     curator (or an assisted lookup step) locates the coordinates on
     the cited assembly/nucleotide record and sets
     `coordinate_provenance: curator_derived`, computing
     `sequence_checksum` at this point. A candidate cannot proceed to
     acceptance with coordinates absent unless
     `coordinate_provenance: not_available` is explicitly set.
3. **Automated validation** (NCBI E-utilities + taxonkit), skipped for
   fields that don't apply to a given `coordinate_provenance`:
   - `accession_resolved` (+ `accession_resolved_version`): resolved
     via the resolver matching the accession's namespace prefix
     (`ncbi_protein:` → NCBI E-utilities; `uniprotkb:` → UniProt REST).
   - `sequence_match`: per gene, percent identity + coverage between
     the segment translation and the cited protein sequence, scored
     pass/warn/fail; record-level `sequence_match.status` is the worst
     of all per-gene statuses.
   - `taxonomy_current`: taxonkit resolves to a current, non-merged
     taxid.
   - Any `fail` or unresolved check keeps `validation.status =
     "needs_review"`. `warn` is surfaced to the reviewer but doesn't
     block review.
4. **Human review gate** (user): review `needs_review` records in
   `db/candidates/` with citations + validation results side by side.
   Accept — move the record directory to
   `db/<Phylum>/<Order_or_Family>/`, set `status="accepted"` (the
   validator enforces this move is atomic with the status change,
   since status and directory location must always agree) — or reject
   — set `status="rejected"` with a required `rejection_reason`, left
   under `db/candidates/`.
5. **GFF3/GBK/protein FASTA generation**: only for accepted records.
6. **Release cut** (periodic, manual trigger): append the current
   `accepted` `record_id@record_version` set to `db/_release.yml`
   under a new release tag, and create a matching git tag in the same
   commit.

Editing an already-`accepted` record (a material change) re-enters at
step 3 in place — `record_version` increments, the record stays under
`db/<Phylum>/`, and `validation.status` drops to `needs_review` until
re-reviewed; it is never moved back to `db/candidates/`.

## Testing / acceptance criteria

- **Schema validation**: every `metadata.yaml` (accepted and
  candidate) validates against `metadata.schema.yaml`; every
  `order.yml` validates against `order.schema.yaml`; idiomorph values
  in `metadata.yaml` are checked against the enum/pattern the matching
  `order.yml` locus entry declares.
- **Cross-file consistency**: for `coordinate_provenance:
  published_explicit`/`curator_derived` records, `locus.gff3`
  coordinates match `metadata.yaml locus.core.segments`; gene
  coordinates/names in `locus.gff3`/`locus.gbk`/`proteins.faa` match
  `genes[]` entries with `present: true`.
- **Status/location invariant**: no record has `validation.status ==
  "accepted"` outside `db/<Phylum>/...`, and none inside it has any
  other status.
- **Edge-case coverage**: per phylum, the seed set includes (where
  literature allows) one `not_available`-coordinate record, one
  `insdc_nucleotide`-sourced record, one record with `present: false`
  genes, and for Basidiomycota, one strain with both an `HD` and a
  `PR` record — each passing its applicable checks.
- **DuckDB build**: build script loads all records without error; row
  counts match file counts; no orphaned foreign keys.
- **Validation script unit tests**: mock NCBI/UniProt responses to
  test accession-resolved / per-gene sequence-match (pass/warn/fail) /
  taxonomy-current logic independently of live network calls.
- **CLI smoke test**: `matpredict curate-db --help` and at least one
  real subcommand run end-to-end against the seed candidate set.
- **Acceptance target**: 5-10 tier-1 records per phylum, including the
  edge cases above, all passing the checks, committed under
  `db/<Phylum>/`, plus working `order.yml`/`order.schema.yaml` for all
  three phyla and a functioning `matpredict curate-db` entry point —
  all four are required for "done," not just the record count.

## Out of scope for this sub-project

- Synteny-based boundary walking (`extended_flank` population beyond
  the schema/empty-list default) — deferred to sub-project 2.
- Tier-2 (homology-only) candidate acceptance.
- Model training / `model_provenance` population — sub-project 2/4.
- Web resource / public browsing UI — sub-project 5.
- `testset/get_MAT.py` unfiltered scraper — superseded by the
  literature-mining workflow above; not modified in this sub-project.
- Full removal of `annotated.yml`/legacy `order.yml` — removed only
  once fully migrated, as a follow-up cleanup.
- `matpredict detect`/`predict`/`classify` subcommands themselves —
  only their shared CLI scaffold (dispatch, config, entry point)
  is built now, as a prerequisite for `curate-db`.
