# MATPredict Sub-project 1: Curated MAT Locus Reference Database

Date: 2026-09-16
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
and `order.yml` (expected gene complement per mating type) stubs.
`testset/Zygo/` has a working prototype for one clade: homology search
(diamond) against known MAT proteins, followed by alignment/tree
building to validate orthology. `testset/get_MAT.py` is an unfiltered
NCBI keyword scraper (searches "mating-type", downloads without
curation) — not part of this design; it predates the curated-DB
approach and is superseded by the workflow below.

## Scope

Populate a curated MAT locus reference database across all three
phyla in parallel (Ascomycota, Basidiomycota, Mucoromycota), restricted
initially to **tier-1 evidence**: loci from peer-reviewed publications
with experimental validation of the locus (e.g. targeted sequencing,
genetic crosses, RFLP), not computational inference alone. Tier-2
(homology-inferred) candidates are tracked in the schema but not
admitted to the accepted database in this phase.

**Done criteria for this sub-project**: schema, curation tooling, and
validation pipeline built and proven on 5-10 tier-1 records per phylum
(~15-30 records total), each passing all validation checks below.
Further population continues afterward as ongoing/parallel work once
sub-project 2 (detection pipeline) starts.

## Directory layout

```
db/
  <Phylum>/
    <Order_or_Family>/
      <taxid>_<genus_species>_<strain>_<idiomorph>/
        locus.gff3          # core MAT locus + flanking genes, on the source assembly's coordinates
        locus.gbk           # full annotated GenBank record for the locus region
        metadata.yaml        # provenance, evidence, taxonomy, validation, versioning
    _schema/
      metadata.schema.yaml   # schema metadata.yaml must satisfy
      duckdb_schema.sql      # DDL for the derived query database
```

Record key: `taxid + strain/isolate + idiomorph`. A heterothallic
species with two idiomorphs gets two records; a homothallic species
with both idiomorphs present gets one record spanning both, per how
the source publication describes the locus.

`Order_or_Family` groups records for browsability as the collection
grows (matches how the eventual web resource would be browsed).

## `metadata.yaml` schema

```yaml
record_id: <taxid>_<strain_slug>_<idiomorph>     # matches directory name
record_version: 1              # increments on any edit to core/extended_flank/genes/evidence
db_release: "2026.09.0"        # dataset-wide release tag, stamped at release-cut time across all records
taxonomy:
  taxid: 4837
  lineage: "k__Fungi;p__Mucoromycota;...;s__Phycomyces_blakesleeanus"  # taxonkit-resolved, cached
  lineage_resolved_date: 2026-09-16
organism:
  species: "Phycomyces blakesleeanus"
  strain: "NRRL 1555"
mating_type:
  idiomorph: "Plus"          # controlled vocabulary per clade (Plus/Minus, MAT1-1/MAT1-2, a/alpha, etc.)
  system: "heterothallic"    # heterothallic | homothallic | pseudohomothallic
locus:
  core:                       # boundary as defined by the source publication; never overwritten in place
    assembly_accession: "GCA_..."
    seq_region: "scaffold_3"
    start: 120345
    end: 128900
    strand: "+"
    definition_note: "core boundary per Idnurm et al. 2008, between tptA and rnhA"
  extended_flank: null        # populated later when synteny work walks boundaries outward; additive, versioned
genes:
  - name: sexP
    protein_accession: "A0A078N0N5_9FUNG"
    role: "core_MAT"          # core_MAT | flanking_conserved | flanking_variable
  - name: tptA
    protein_accession: "B0F2G7_PHYBL"
    role: "flanking_conserved"
evidence:
  tier: 1                     # 1 = published + experimentally validated (only tier admitted to v1 DB)
  citation:
    pmid: "18248337"
    doi: "10.1128/EC.00281-08"
  experimental_method: "targeted sequencing + genetic crosses"
validation:
  accession_resolved: true
  accession_resolved_date: 2026-09-16
  sequence_match: true         # translated CDS matches cited protein accession
  sequence_match_notes: ""
  taxonomy_current: true
  status: "accepted"           # accepted | needs_review | rejected
curation:
  proposed_by: "literature-mining-agent"
  reviewed_by: "jason.stajich@ucr.edu"
  reviewed_date: 2026-09-16
model_provenance: null         # reserved: populated by sub-project 2/4 when a model is trained citing this record
```

Fields not yet known at curation time (e.g. an accession the paper
doesn't specify) are left blank, never inferred or guessed.

## DuckDB schema (query cache, not source of truth)

Generated by a build script that walks all `metadata.yaml` files.
GFF3/GBK/YAML on disk remain authoritative; DuckDB is rebuilt from them,
never edited directly.

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
    strain TEXT,
    idiomorph TEXT NOT NULL,
    mating_system TEXT,                -- heterothallic | homothallic | pseudohomothallic
    phylum TEXT NOT NULL,
    order_or_family TEXT,
    record_version INTEGER NOT NULL,
    db_release TEXT NOT NULL,
    assembly_accession TEXT,
    seq_region TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    strand TEXT,
    definition_note TEXT,
    evidence_tier INTEGER NOT NULL,
    pmid TEXT,
    doi TEXT,
    experimental_method TEXT,
    validation_status TEXT NOT NULL,   -- accepted | needs_review | rejected
    accession_resolved BOOLEAN,
    sequence_match BOOLEAN,
    taxonomy_current BOOLEAN,
    proposed_by TEXT,
    reviewed_by TEXT,
    reviewed_date DATE,
    gff3_path TEXT NOT NULL,
    gbk_path TEXT NOT NULL,
    metadata_path TEXT NOT NULL
);

CREATE TABLE locus_gene (
    record_id TEXT REFERENCES locus_record(record_id),
    gene_name TEXT NOT NULL,
    protein_accession TEXT,
    role TEXT NOT NULL,                -- core_MAT | flanking_conserved | flanking_variable
    PRIMARY KEY (record_id, gene_name)
);

CREATE TABLE extended_flank (
    record_id TEXT REFERENCES locus_record(record_id),
    walked_out_date DATE,
    seq_region TEXT,
    start_pos INTEGER,
    end_pos INTEGER,
    method TEXT,                       -- how the extension was derived (synteny tool, manual, etc.)
    notes TEXT
);
```

## Curation workflow

Enforces verifiability at every step; no record reaches `accepted`
without passing automated validation and explicit human sign-off.

1. **Literature search** (per phylum/order, via PubMed): agent
   proposes candidates. Every candidate must carry a resolvable
   PMID/DOI and the specific text/table the claim was drawn from. No
   candidate is proposed without one.
2. **Candidate extraction** → draft `metadata.yaml` with
   `evidence.tier=1`, `validation.status="needs_review"`. Fields not
   explicitly stated in the source are left blank, never inferred.
3. **Automated validation** (NCBI E-utilities + taxonkit):
   - `accession_resolved`: does the cited accession still resolve
     (not suppressed/replaced)?
   - `sequence_match`: does the assembly region translate to / contain
     the cited protein sequence?
   - `taxonomy_current`: does taxonkit resolve the species to a
     current, non-merged taxid?
   - Any failure keeps `validation.status = "needs_review"`; no
     auto-accept.
4. **Human review gate** (user): review `needs_review` records with
   citation + validation results side by side. Accept
   (`status="accepted"`, enters the current `db_release`) or reject
   (`status="rejected"` with a reason, kept for traceability so the
   mining agent doesn't re-propose it blindly).
5. **GFF3/GBK generation**: only for accepted records, built from the
   validated coordinates/accessions.

## Testing / acceptance criteria

- **Schema validation**: every `metadata.yaml` validates against
  `metadata.schema.yaml` (required fields present, controlled
  vocabulary enforced for `idiomorph`, `role`, `validation.status`).
- **Cross-file consistency**: `locus.gff3` coordinates match
  `metadata.yaml locus.core`; gene names in `locus.gff3`/`locus.gbk`
  match `genes[].name`.
- **DuckDB build**: build script loads all records without error; row
  counts match file counts.
- **Validation script unit tests**: mock NCBI responses to test
  accession-resolved / sequence-match / taxonomy-current logic
  independently of live network calls.
- **Acceptance target**: 5-10 tier-1 records per phylum, all passing
  the above, committed to `db/`.

## Out of scope for this sub-project

- Synteny-based boundary walking (`extended_flank` population) —
  schema supports it, but population is deferred to sub-project 2.
- Tier-2 (homology-only) candidate acceptance.
- Model training / `model_provenance` population — sub-project 2/4.
- Web resource / public browsing UI — sub-project 5.
- `testset/get_MAT.py` unfiltered scraper — superseded by the
  literature-mining workflow above; not modified in this sub-project.
