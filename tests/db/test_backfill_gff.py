from __future__ import annotations

from unittest.mock import MagicMock

import yaml

from MATPredict.db.cli import build_gff_for_record, find_records_missing_proteins_faa


def _write_record(db_root, phylum, order_or_family, record_id, genes=None):
    record_dir = db_root / phylum / order_or_family / record_id
    record_dir.mkdir(parents=True)
    metadata = {
        "record_id": record_id,
        "genes": genes if genes is not None else [
            {"gene_index": 0, "name": "G1", "present": True, "protein_accession": "ncbi_protein:ABC1.1"},
        ],
    }
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(metadata))
    return record_dir


def test_find_records_missing_proteins_faa_skips_records_that_already_have_one(tmp_path):
    db_root = tmp_path / "db"
    _write_record(db_root, "Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")
    complete_dir = _write_record(db_root, "Ascomycota", "Eurotiales", "222_b_MAT_MAT1-2")
    (complete_dir / "proteins.faa").write_text(">already|gene_index=0|name=G1|role=core_MAT\nMSEQ\n")

    missing = find_records_missing_proteins_faa(db_root)

    assert missing == [("Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")]


def test_find_records_missing_proteins_faa_excludes_candidates(tmp_path):
    db_root = tmp_path / "db"
    _write_record(db_root, "candidates", "Ascomycota", "999_c_MAT_MAT1-1")

    assert find_records_missing_proteins_faa(db_root) == []


def test_backfill_missing_proteins_faa_isolates_one_record_failure(tmp_path, monkeypatch):
    db_root = tmp_path / "db"
    _write_record(db_root, "Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")
    _write_record(db_root, "Ascomycota", "Onygenales", "222_b_MAT_MAT1-2")

    calls = []

    def fake_build(db_root, phylum, order_or_family, record_id, ncbi, uniprot):
        calls.append(record_id)
        if record_id == "111_a_MAT_MAT1-1":
            raise ConnectionError("simulated NCBI outage")

    monkeypatch.setattr("MATPredict.db.cli.build_gff_for_record", fake_build)

    from MATPredict.db.cli import backfill_missing_proteins_faa

    succeeded, failed = backfill_missing_proteins_faa(db_root, ncbi=MagicMock(), uniprot=MagicMock())

    assert calls == ["111_a_MAT_MAT1-1", "222_b_MAT_MAT1-2"]  # second record still attempted
    assert succeeded == [("Ascomycota", "Onygenales", "222_b_MAT_MAT1-2")]
    assert len(failed) == 1
    assert failed[0][0] == ("Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")
    assert "simulated NCBI outage" in failed[0][1]
