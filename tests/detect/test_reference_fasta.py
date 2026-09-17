"""Test reference protein FASTA builder."""
from __future__ import annotations

from MATPredict.detect.reference_fasta import build_reference_fasta


def test_build_reference_fasta_concatenates_and_rewrites_headers(tmp_path):
    """Test that build_reference_fasta transforms gff_export headers correctly."""
    record_dir = tmp_path / "Basidiomycota" / "Ustilaginaceae" / "5270_521_aLocus_a1"
    record_dir.mkdir(parents=True)
    (record_dir / "proteins.faa").write_text(
        ">5270_521_aLocus_a1|gene_index=0|name=mfa1|role=core_MAT\nMKV\n"
    )
    out = build_reference_fasta(tmp_path, tmp_path / "combined.faa")
    text = out.read_text()
    assert text.startswith(">5270_521_aLocus_a1|gene0|mfa1\n")
    assert "MKV" in text


def test_build_reference_fasta_excludes_candidates(tmp_path):
    """Finding E regression: db/candidates/ holds not-yet-accepted (needs_review)
    records, which family_registry.load_record_families and
    pipeline._short_orf_genes both already exclude as non-authoritative.
    Including a candidate's proteins here would let it be hit by search, only
    for search._attribute to silently drop the hit later (its record_id is
    absent from the accepted-records index) -- wasted search cost with no
    signal to the user."""
    accepted_dir = tmp_path / "Basidiomycota" / "Ustilaginaceae" / "accepted_record"
    accepted_dir.mkdir(parents=True)
    (accepted_dir / "proteins.faa").write_text(
        ">accepted_record|gene_index=0|name=mfa1|role=core_MAT\nMKV\n"
    )
    candidate_dir = tmp_path / "candidates" / "Ascomycota" / "candidate_record"
    candidate_dir.mkdir(parents=True)
    (candidate_dir / "proteins.faa").write_text(
        ">candidate_record|gene_index=0|name=mfa1|role=core_MAT\nQQQ\n"
    )
    out = build_reference_fasta(tmp_path, tmp_path / "combined.faa")
    text = out.read_text()
    assert "accepted_record" in text
    assert "candidate_record" not in text
    assert "QQQ" not in text


def test_build_reference_fasta_logs_malformed_headers(tmp_path, caplog):
    """Test that malformed headers are skipped and a warning is logged."""
    record_dir = tmp_path / "Basidiomycota" / "Ustilaginaceae" / "test_record"
    record_dir.mkdir(parents=True)
    faa_file = record_dir / "proteins.faa"
    faa_file.write_text(
        ">test_record|gene_index=0|name=good_gene|role=core_MAT\nMKV\n"
        ">malformed_header\nFAA\n"
    )
    out = build_reference_fasta(tmp_path, tmp_path / "combined.faa")
    text = out.read_text()
    # Well-formed header should be included
    assert ">test_record|gene0|good_gene" in text
    assert "MKV" in text
    # Malformed header should not be included
    assert "malformed_header" not in text
    assert "FAA" not in text
    # Warning should be logged
    assert "Skipping malformed header" in caplog.text
    assert "malformed_header" in caplog.text
