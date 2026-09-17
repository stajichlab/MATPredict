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
