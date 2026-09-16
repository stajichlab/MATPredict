from __future__ import annotations

import copy
import shutil
from pathlib import Path

import pytest
import yaml

from MATPredict.db.curate import CurationError, accept_candidate, propose_candidate, reject_candidate

_REAL_MUCOROMYCOTA_ORDER_YML = Path(__file__).resolve().parents[2] / "db" / "Mucoromycota" / "order.yml"


@pytest.fixture(autouse=True)
def _seed_order_yml(tmp_path):
    """propose_candidate/accept_candidate validate mating_type.idiomorphs against
    db/<phylum>/order.yml, so every test needs a copy of the real Mucoromycota order.yml
    under its tmp_path db_root."""
    dest_dir = tmp_path / "Mucoromycota"
    dest_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(_REAL_MUCOROMYCOTA_ORDER_YML, dest_dir / "order.yml")


RECORD = {
    "record_id": "4837_nrrl-1555_MAT_Plus",
    "record_version": 1,
    "taxonomy": {"taxid": 4837, "lineage": "k__Fungi;p__Mucoromycota", "lineage_resolved_date": "2026-09-16"},
    "organism": {"species": "Phycomyces blakesleeanus", "strain": {"name": "NRRL 1555", "known": True}},
    "mating_type": {"locus_name": "MAT", "idiomorphs": ["Plus"], "system": "heterothallic"},
    "locus": {"coordinate_provenance": "not_available", "excluded_from_coordinate_benchmark": True},
    "genes": [{"gene_index": 0, "name": "sexP", "protein_accession": "ncbi_protein:AAB12345.1", "role": "core_MAT", "present": True}],
    "evidence": {
        "locus_existence": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "boundaries": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "idiomorph_assignment": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
    },
    "validation": {"status": "needs_review", "rejection_reason": None},
    "curation": {"proposed_by": "literature-mining-agent", "proposal_dedupe_key": "18248337|4837|MAT"},
}


def test_propose_writes_candidate_dir(tmp_path):
    candidate_dir = propose_candidate(tmp_path, "Mucoromycota", RECORD)
    assert candidate_dir == tmp_path / "candidates" / "Mucoromycota" / RECORD["record_id"]
    written = yaml.safe_load((candidate_dir / "metadata.yaml").read_text())
    assert written["validation"]["status"] == "needs_review"


def test_propose_rejects_invalid_record(tmp_path):
    bad = copy.deepcopy(RECORD)
    del bad["evidence"]
    with pytest.raises(CurationError):
        propose_candidate(tmp_path, "Mucoromycota", bad)


def test_accept_moves_directory_and_sets_status(tmp_path):
    propose_candidate(tmp_path, "Mucoromycota", RECORD)
    accepted_dir = accept_candidate(tmp_path, "Mucoromycota", "Mucorales", RECORD["record_id"])
    assert accepted_dir == tmp_path / "Mucoromycota" / "Mucorales" / RECORD["record_id"]
    assert not (tmp_path / "candidates" / "Mucoromycota" / RECORD["record_id"]).exists()
    written = yaml.safe_load((accepted_dir / "metadata.yaml").read_text())
    assert written["validation"]["status"] == "accepted"


def test_reject_sets_status_and_reason_in_place(tmp_path):
    propose_candidate(tmp_path, "Mucoromycota", RECORD)
    reject_candidate(tmp_path, "Mucoromycota", RECORD["record_id"], reason="accession no longer resolves")
    written_path = tmp_path / "candidates" / "Mucoromycota" / RECORD["record_id"] / "metadata.yaml"
    written = yaml.safe_load(written_path.read_text())
    assert written["validation"]["status"] == "rejected"
    assert written["validation"]["rejection_reason"] == "accession no longer resolves"


def test_reject_requires_a_reason(tmp_path):
    propose_candidate(tmp_path, "Mucoromycota", RECORD)
    with pytest.raises(CurationError):
        reject_candidate(tmp_path, "Mucoromycota", RECORD["record_id"], reason="")


# --- Ruling 1: schema.validate_idiomorphs wired into propose/accept ---


def test_propose_rejects_idiomorph_not_in_order_vocabulary(tmp_path):
    bad = copy.deepcopy(RECORD)
    bad["mating_type"] = {"locus_name": "MAT", "idiomorphs": ["NotAValue"], "system": "heterothallic"}
    with pytest.raises(CurationError):
        propose_candidate(tmp_path, "Mucoromycota", bad)


def test_propose_raises_curation_error_not_keyerror_when_mating_type_missing(tmp_path):
    bad = copy.deepcopy(RECORD)
    del bad["mating_type"]
    with pytest.raises(CurationError):
        propose_candidate(tmp_path, "Mucoromycota", bad)


# --- Ruling 2: proposal_dedupe_key duplicate-against-rejected check ---


def test_propose_rejects_duplicate_dedupe_key_against_rejected_candidate(tmp_path):
    propose_candidate(tmp_path, "Mucoromycota", RECORD)
    reject_candidate(tmp_path, "Mucoromycota", RECORD["record_id"], reason="accession no longer resolves")

    new_record = copy.deepcopy(RECORD)
    new_record["record_id"] = "4837_nrrl-1555_MAT_Minus"
    new_record["mating_type"] = {"locus_name": "MAT", "idiomorphs": ["Minus"], "system": "heterothallic"}
    # same proposal_dedupe_key as the rejected candidate above
    new_record["curation"]["proposal_dedupe_key"] = "18248337|4837|MAT"

    with pytest.raises(CurationError):
        propose_candidate(tmp_path, "Mucoromycota", new_record)
