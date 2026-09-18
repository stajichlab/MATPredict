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
