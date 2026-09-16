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
