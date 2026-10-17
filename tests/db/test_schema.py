from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from MATPredict.db import schema

VALID_ORDER = {
    "phylum": "Mucoromycota",
    "loci": [
        {
            "locus_name": "MAT",
            "vocabulary_type": "enum",
            "idiomorph_values": ["Plus", "Minus"],
            "taxonomic_scope": [4761],
            "genes": [
                {"name": "tptA", "role": "flanking_conserved"},
                {"name": "sexP", "role": "core_MAT", "present_in_idiomorphs": ["Plus"]},
                {"name": "sexM", "role": "core_MAT", "present_in_idiomorphs": ["Minus"]},
                {"name": "rnhA", "role": "flanking_conserved"},
            ],
        }
    ],
}

VALID_RECORD = {
    "record_id": "4837_nrrl-1555_MAT_Plus",
    "record_version": 1,
    "taxonomy": {"taxid": 4837, "lineage": "k__Fungi;p__Mucoromycota", "lineage_resolved_date": "2026-09-16"},
    "organism": {"species": "Phycomyces blakesleeanus", "strain": {"name": "NRRL 1555", "known": True}},
    "mating_type": {"locus_name": "MAT", "idiomorphs": ["Plus"], "system": "heterothallic"},
    "locus": {"coordinate_provenance": "not_available", "excluded_from_coordinate_benchmark": True},
    "genes": [
        {"gene_index": 0, "name": "sexP", "protein_accession": "ncbi_protein:AAB12345.1",
         "role": "core_MAT", "present": True},
    ],
    "evidence": {
        "locus_existence": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "boundaries": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "idiomorph_assignment": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
    },
    "validation": {"status": "needs_review", "rejection_reason": None},
    "curation": {"proposed_by": "literature-mining-agent"},
}


def test_valid_order_document_passes():
    assert schema.validate_order(VALID_ORDER) == []


def test_valid_record_passes():
    assert schema.validate_metadata(VALID_RECORD) == []


def test_record_missing_required_field_fails():
    bad = copy.deepcopy(VALID_RECORD)
    del bad["evidence"]
    errors = schema.validate_metadata(bad)
    assert errors
    assert any("evidence" in e for e in errors)


def test_idiomorph_enum_rejects_unknown_value():
    bad = copy.deepcopy(VALID_RECORD)
    bad["mating_type"]["idiomorphs"] = ["NotAValue"]
    errors = schema.validate_idiomorphs(bad, VALID_ORDER)
    assert errors


def test_idiomorph_pattern_locus_accepts_allele_ids():
    pattern_order = {
        "phylum": "Basidiomycota",
        "loci": [{"locus_name": "HD", "vocabulary_type": "pattern", "idiomorph_pattern": "^A[0-9]+$",
                   "taxonomic_scope": [4982],
                   "genes": [{"name": "HD1", "role": "core_MAT"}]}],
    }
    record = copy.deepcopy(VALID_RECORD)
    record["mating_type"] = {"locus_name": "HD", "idiomorphs": ["A1"], "system": "heterothallic"}
    assert schema.validate_idiomorphs(record, pattern_order) == []


def test_validate_order_requires_taxonomic_scope():
    doc = {
        "phylum": "TestPhylum",
        "loci": [
            {
                "locus_name": "MAT",
                "vocabulary_type": "enum",
                "idiomorph_values": ["a", "alpha"],
                "genes": [{"name": "STE3", "role": "core_MAT"}],
            }
        ],
    }
    errors = schema.validate_order(doc)
    assert any("taxonomic_scope" in e for e in errors)


def test_validate_order_accepts_taxonomic_scope():
    doc = {
        "phylum": "TestPhylum",
        "loci": [
            {
                "locus_name": "MAT",
                "vocabulary_type": "enum",
                "idiomorph_values": ["a", "alpha"],
                "taxonomic_scope": [4930],
                "genes": [{"name": "STE3", "role": "core_MAT"}],
            }
        ],
    }
    assert schema.validate_order(doc) == []


def test_record_valid_without_exon_fields():
    """Pre-existing gene shape (no exons/codon_start/transl_table) must still validate."""
    assert schema.validate_metadata(VALID_RECORD) == []


def test_record_valid_with_exon_fields():
    record = copy.deepcopy(VALID_RECORD)
    record["genes"][0]["exons"] = [{"start": 100, "end": 150}, {"start": 200, "end": 260}]
    record["genes"][0]["codon_start"] = 1
    record["genes"][0]["transl_table"] = 1
    assert schema.validate_metadata(record) == []


@pytest.mark.parametrize("bad_codon_start", [0, 4])
def test_record_rejects_invalid_codon_start(bad_codon_start):
    record = copy.deepcopy(VALID_RECORD)
    record["genes"][0]["codon_start"] = bad_codon_start
    errors = schema.validate_metadata(record)
    assert errors
    assert any("codon_start" in e for e in errors)


@pytest.mark.parametrize("phylum", ["Ascomycota", "Basidiomycota", "Mucoromycota"])
def test_real_order_yml_files_validate(phylum):
    repo_root = Path(__file__).resolve().parents[2]
    order_path = repo_root / "db" / phylum / "order.yml"
    doc = yaml.safe_load(order_path.read_text())
    assert schema.validate_order(doc) == []


# --- gene-vocabulary reconciliation (order.yml `genes` vs curated record `genes`) ---


def test_gene_vocabulary_accepts_a_declared_name():
    assert schema.validate_gene_vocabulary(VALID_RECORD, VALID_ORDER) == []


def test_gene_vocabulary_rejects_an_undeclared_name():
    bad = copy.deepcopy(VALID_RECORD)
    bad["genes"][0]["name"] = "sexP_variant7"
    errors = schema.validate_gene_vocabulary(bad, VALID_ORDER)
    assert len(errors) == 1
    assert "sexP_variant7" in errors[0]


def test_gene_vocabulary_skips_absent_genes():
    absent = copy.deepcopy(VALID_RECORD)
    absent["genes"][0]["name"] = "sexP_variant7"
    absent["genes"][0]["present"] = False
    assert schema.validate_gene_vocabulary(absent, VALID_ORDER) == []


def test_gene_vocabulary_reports_an_unknown_locus():
    unknown = copy.deepcopy(VALID_RECORD)
    unknown["mating_type"]["locus_name"] = "NotALocus"
    assert schema.validate_gene_vocabulary(unknown, VALID_ORDER) == ["no order.yml locus entry named 'NotALocus'"]


def _accepted_records():
    """Every accepted curated record, paired with its phylum's order.yml document.

    `db/candidates/` is excluded the same way `family_registry.load_record_families`
    excludes it: those are proposed, not accepted, records.
    """
    repo_root = Path(__file__).resolve().parents[2]
    db_root = repo_root / "db"
    orders: dict[str, dict] = {}
    for meta_path in sorted(db_root.glob("*/*/*/metadata.yaml")):
        phylum = meta_path.relative_to(db_root).parts[0]
        if phylum == "candidates":
            continue
        if phylum not in orders:
            orders[phylum] = yaml.safe_load((db_root / phylum / "order.yml").read_text())
        yield meta_path, orders[phylum]


def test_every_accepted_record_gene_name_is_declared_by_its_family():
    """Regression guard for the silent-orphaning bug.

    `detect.search._attribute` drops any hit whose curated gene name is not in its
    family's `order.yml` `genes` list. Before this guard, 15 of 95 accepted curated
    reference proteins were unreachable that way (Basidiomycota PR, Balpha and Bbeta),
    which made Balpha and Bbeta undetectable in every genome including their own.
    """
    failures = []
    for meta_path, order_doc in _accepted_records():
        record = yaml.safe_load(meta_path.read_text())
        for error in schema.validate_gene_vocabulary(record, order_doc):
            failures.append(f"{meta_path}: {error}")
    assert failures == []


def test_every_gene_class_used_in_order_yml_is_a_member_of_the_schema_enum():
    """Permanent net against a typo'd gene_class silently disappearing from figures.

    `write_genbank` joins a record's genes to their phylum's `order.yml` to emit the
    `/gene_class` CDS qualifier, and `draw.py`/`synteny.py` colour and label by it.
    A misspelled class does not raise anywhere -- it simply produces a class nothing
    else shares, so the gene quietly gets its own colour and its own label. The enum
    is read from `db/_schema/order.schema.yaml` itself rather than restated here, so
    this test cannot drift from the schema it guards.
    """
    allowed = set(
        schema.load_order_schema()["properties"]["loci"]["items"]["properties"]["genes"]
        ["items"]["properties"]["gene_class"]["enum"]
    )
    repo_root = Path(__file__).resolve().parents[2]
    offenders = []
    for order_path in sorted((repo_root / "db").glob("*/order.yml")):
        doc = yaml.safe_load(order_path.read_text()) or {}
        for locus in doc.get("loci", []):
            for gene in locus.get("genes", []):
                gene_class = gene.get("gene_class")
                if gene_class is not None and gene_class not in allowed:
                    offenders.append(
                        f"{order_path.name}:{locus['locus_name']}:{gene['name']} -> {gene_class}"
                    )
    assert offenders == [], f"gene_class values outside the schema enum: {offenders}"


@pytest.mark.parametrize("system", ["heterothallic", "homothallic", "pseudohomothallic", "unisexual_switching"])
def test_mating_system_enum_accepts_known_values(system):
    record = copy.deepcopy(VALID_RECORD)
    record["mating_type"]["system"] = system
    assert schema.validate_metadata(record) == []


def test_mating_system_enum_rejects_unknown_value():
    record = copy.deepcopy(VALID_RECORD)
    record["mating_type"]["system"] = "unisexual"
    errors = schema.validate_metadata(record)
    assert errors and any("system" in e for e in errors)
# --- gene-level `status` (remnant|pseudogene) ---------------------------------


def _with_gene(**fields):
    record = copy.deepcopy(VALID_RECORD)
    record["genes"][0].update(fields)
    return record


@pytest.mark.parametrize("status", ["remnant", "pseudogene"])
def test_schema_accepts_gene_status(status):
    assert schema.validate_metadata(_with_gene(status=status, remnant_note="x")) == []


@pytest.mark.parametrize("status", ["partial", "intact", "Remnant", "", None])
def test_schema_rejects_other_gene_status_values_and_partial(status):
    assert schema.validate_metadata(_with_gene(status=status)) != []


def test_status_is_optional_and_absent_means_intact():
    assert "status" not in VALID_RECORD["genes"][0]
    assert schema.validate_metadata(VALID_RECORD) == []
    assert schema.validate_gene_status(VALID_RECORD) == []


def test_validate_gene_status_present_false_with_no_model_is_valid():
    gene = dict(protein_accession=None, present=False, status="remnant", start=None, end=None)
    assert schema.validate_gene_status(_with_gene(**gene)) == []


def test_validate_gene_status_present_false_with_a_model_is_flagged():
    errors = schema.validate_gene_status(_with_gene(present=False, status="remnant"))
    assert len(errors) == 1 and "present: false" in errors[0]


def test_validate_gene_status_present_true_with_a_model_is_valid():
    assert schema.validate_gene_status(_with_gene(status="pseudogene")) == []


def test_validate_gene_status_present_true_without_a_model_is_flagged():
    errors = schema.validate_gene_status(_with_gene(status="remnant", protein_accession=None))
    assert len(errors) == 1 and "records no model" in errors[0]


def test_validate_gene_status_flags_an_orphan_remnant_note():
    errors = schema.validate_gene_status(_with_gene(remnant_note="truncated"))
    assert len(errors) == 1 and "no status" in errors[0]


# Records migrated to carry a gene `status` (decision paper
# docs/notes/2026-10-17_remnant-status.md, option B step 1): (record_id, gene, status, present).
# sxlc146 MAT1-1-2 is deliberately absent: its 391-aa protein is full length against
# Diaporthales references (397, 406 aa), so it is not labelled.
MIGRATED_STATUS_GENES = [
    ("40559_sas405_MAT_MAT1-2", "MAT1-1-1", "remnant", False),
    ("28548_06cwm-g22_MAT_combined", "MAT1-1-1", "remnant", False),
    ("1708542_exf-10342_wallMAT_v2", "HMG", "remnant", False),
    ("746128_a1163_MAT_MAT1-1", "MAT1-2-1", "remnant", True),
    ("230073_uamh-11059_MAT_MAT1-2", "MAT1-1-1", "remnant", True),
    ("107463_unknown-1_MAT_MAT1-2", "MAT1-2-1", "pseudogene", True),
]


def test_migrated_records_carry_the_gene_status_flag():
    records = {}
    for meta_path, _order in _accepted_records():
        record = yaml.safe_load(meta_path.read_text())
        records[record["record_id"]] = record
    for record_id, gene_name, status, present in MIGRATED_STATUS_GENES:
        genes = [g for g in records[record_id]["genes"] if g["name"] == gene_name]
        assert len(genes) == 1, record_id
        assert genes[0].get("status") == status, record_id
        assert genes[0]["present"] is present, record_id
        assert records[record_id]["record_version"] >= 2, record_id
    # exactly these genes carry a status, nothing else in the database
    flagged = {
        (rid, g["name"]) for rid, r in records.items() for g in r["genes"] if g.get("status")
    }
    assert flagged == {(rid, gene) for rid, gene, _s, _p in MIGRATED_STATUS_GENES}


def test_every_accepted_record_validates_against_the_schema_and_gene_status_rules():
    failures = []
    for meta_path, _order in _accepted_records():
        record = yaml.safe_load(meta_path.read_text())
        for error in schema.validate_metadata(record) + schema.validate_gene_status(record):
            failures.append(f"{meta_path}: {error}")
    assert failures == []


def test_no_status_gene_protein_is_in_a_proteins_faa_reference_file():
    """Step 2: proteins.faa is what detect reads as the reference set, so it must not
    carry a protein for a gene whose metadata says remnant|pseudogene."""
    failures = []
    for meta_path, _order in _accepted_records():
        record = yaml.safe_load(meta_path.read_text())
        flagged = {g["gene_index"] for g in record["genes"] if g.get("status")}
        faa = meta_path.parent / "proteins.faa"
        if not flagged or not faa.exists():
            continue
        for line in faa.read_text().splitlines():
            if line.startswith(">"):
                fields = dict(p.split("=", 1) for p in line.split("|")[1:] if "=" in p)
                if int(fields["gene_index"]) in flagged:
                    failures.append(f"{faa}: {line}")
    assert failures == []
