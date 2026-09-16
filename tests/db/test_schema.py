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
                   "genes": [{"name": "HD1", "role": "core_MAT"}]}],
    }
    record = copy.deepcopy(VALID_RECORD)
    record["mating_type"] = {"locus_name": "HD", "idiomorphs": ["A1"], "system": "heterothallic"}
    assert schema.validate_idiomorphs(record, pattern_order) == []


@pytest.mark.parametrize("phylum", ["Ascomycota", "Basidiomycota", "Mucoromycota"])
def test_real_order_yml_files_validate(phylum):
    repo_root = Path(__file__).resolve().parents[2]
    order_path = repo_root / "db" / phylum / "order.yml"
    doc = yaml.safe_load(order_path.read_text())
    assert schema.validate_order(doc) == []
