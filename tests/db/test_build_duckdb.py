from __future__ import annotations

import copy

import duckdb
import yaml

from MATPredict.db.build_duckdb import build

RECORD = {
    "record_id": "4837_nrrl-1555_MAT_Plus",
    "record_version": 1,
    "taxonomy": {"taxid": 4837, "lineage": "k__Fungi;p__Mucoromycota;s__Phycomyces_blakesleeanus", "lineage_resolved_date": "2026-09-16"},
    "organism": {"species": "Phycomyces blakesleeanus", "strain": {"name": "NRRL 1555", "known": True}},
    "mating_type": {"locus_name": "MAT", "idiomorphs": ["Plus"], "system": "heterothallic"},
    "locus": {"coordinate_provenance": "not_available", "excluded_from_coordinate_benchmark": True},
    "genes": [{"gene_index": 0, "name": "sexP", "protein_accession": "ncbi_protein:AAB12345.1", "role": "core_MAT", "present": True}],
    "evidence": {
        "locus_existence": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "boundaries": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "idiomorph_assignment": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
    },
    "validation": {"status": "accepted", "rejection_reason": None},
    "curation": {"proposed_by": "literature-mining-agent"},
}


def test_build_loads_all_records(tmp_path):
    record_dir = tmp_path / "Mucoromycota" / "Mucorales" / RECORD["record_id"]
    record_dir.mkdir(parents=True)
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(RECORD))

    out_path = tmp_path / "matpredict.duckdb"
    build(db_root=tmp_path, out_path=out_path)

    con = duckdb.connect(str(out_path))
    count = con.execute("SELECT COUNT(*) FROM locus_record").fetchone()[0]
    assert count == 1
    row = con.execute("SELECT phylum, validation_status FROM locus_record").fetchone()
    assert row == ("Mucoromycota", "accepted")
    gene_count = con.execute("SELECT COUNT(*) FROM locus_gene").fetchone()[0]
    assert gene_count == 1
    con.close()


def test_accepted_record_with_generated_files_gets_real_paths(tmp_path):
    record = copy.deepcopy(RECORD)
    record["record_id"] = "4837_nrrl-1555_MAT_Plus_files"
    record_dir = tmp_path / "Mucoromycota" / "Mucorales" / record["record_id"]
    record_dir.mkdir(parents=True)
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(record))
    (record_dir / "locus.gff3").write_text("")
    (record_dir / "locus.gbk").write_text("")
    (record_dir / "proteins.faa").write_text("")

    out_path = tmp_path / "matpredict.duckdb"
    build(db_root=tmp_path, out_path=out_path)

    con = duckdb.connect(str(out_path))
    row = con.execute(
        "SELECT gff3_path, gbk_path, proteins_fasta_path FROM locus_record WHERE record_id = ?",
        [record["record_id"]],
    ).fetchone()
    con.close()

    assert row[0] == str(record_dir / "locus.gff3")
    assert row[1] == str(record_dir / "locus.gbk")
    assert row[2] == str(record_dir / "proteins.faa")


def test_accepted_record_without_generated_files_gets_null_paths(tmp_path):
    record = copy.deepcopy(RECORD)
    record["record_id"] = "4837_nrrl-1555_MAT_Plus_nofiles"
    record_dir = tmp_path / "Mucoromycota" / "Mucorales" / record["record_id"]
    record_dir.mkdir(parents=True)
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(record))
    # Deliberately do not create locus.gff3 / locus.gbk / proteins.faa.

    out_path = tmp_path / "matpredict.duckdb"
    build(db_root=tmp_path, out_path=out_path)

    con = duckdb.connect(str(out_path))
    row = con.execute(
        "SELECT gff3_path, gbk_path, proteins_fasta_path FROM locus_record WHERE record_id = ?",
        [record["record_id"]],
    ).fetchone()
    con.close()

    assert row == (None, None, None)


def test_locus_idiomorph_gets_one_row_per_idiomorph_value(tmp_path):
    record = copy.deepcopy(RECORD)
    record["record_id"] = "4837_nrrl-1555_MAT_combined"
    record["mating_type"] = {"locus_name": "MAT", "idiomorphs": ["Plus", "Minus"], "system": "homothallic"}
    record_dir = tmp_path / "Mucoromycota" / "Mucorales" / record["record_id"]
    record_dir.mkdir(parents=True)
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(record))

    out_path = tmp_path / "matpredict.duckdb"
    build(db_root=tmp_path, out_path=out_path)

    con = duckdb.connect(str(out_path))
    rows = con.execute(
        "SELECT idiomorph_value FROM locus_idiomorph WHERE record_id = ? ORDER BY idiomorph_value",
        [record["record_id"]],
    ).fetchall()
    con.close()

    assert rows == [("Minus",), ("Plus",)]


def test_candidate_record_gets_null_paths_even_with_files_present(tmp_path):
    record = copy.deepcopy(RECORD)
    record["record_id"] = "4837_nrrl-1555_MAT_Plus_candidate"
    record["validation"] = {"status": "needs_review", "rejection_reason": None}
    record_dir = tmp_path / "candidates" / "Mucoromycota" / "Mucorales" / record["record_id"]
    record_dir.mkdir(parents=True)
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(record))
    # A candidate should never have these generated files, but even if
    # they happen to be present on disk, a candidate must never get
    # non-null paths.
    (record_dir / "locus.gff3").write_text("")
    (record_dir / "locus.gbk").write_text("")
    (record_dir / "proteins.faa").write_text("")

    out_path = tmp_path / "matpredict.duckdb"
    build(db_root=tmp_path, out_path=out_path)

    con = duckdb.connect(str(out_path))
    row = con.execute(
        "SELECT gff3_path, gbk_path, proteins_fasta_path FROM locus_record WHERE record_id = ?",
        [record["record_id"]],
    ).fetchone()
    con.close()

    assert row == (None, None, None)
