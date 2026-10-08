"""`matpredict detect` says what to do when the genome file is compressed (seam: the CLI entry point)."""
from __future__ import annotations

import gzip
import logging

import pytest

from MATPredict.__main__ import main


def _run(genome, out):
    return main(["detect", "--genome", str(genome), "--out-dir", str(out), "--taxid", "5507"])


def test_gzipped_genome_is_refused_with_a_decompress_hint(tmp_path, caplog):
    g = tmp_path / "asm.fa.gz"
    with gzip.open(g, "wt") as fh:
        fh.write(">c1\nACGT\n")
    with caplog.at_level(logging.ERROR):
        assert _run(g, tmp_path / "out") == 1
    msg = " ".join(r.getMessage() for r in caplog.records)
    assert "asm.fa.gz" in msg and "gzip" in msg.lower() and "decompress" in msg.lower()
    assert not (tmp_path / "out").exists()          # refused before any output was created


def test_compression_is_recognised_by_content_not_by_name(tmp_path, caplog):
    g = tmp_path / "asm.fna"                         # gzip data under a plain name
    with gzip.open(g, "wt") as fh:
        fh.write(">c1\nACGT\n")
    with caplog.at_level(logging.ERROR):
        assert _run(g, tmp_path / "out") == 1
    assert "decompress" in " ".join(r.getMessage() for r in caplog.records).lower()


def test_zstd_genome_is_refused_with_a_decompress_hint(tmp_path, caplog):
    g = tmp_path / "asm.fa.zst"
    g.write_bytes(b"\x28\xb5\x2f\xfd" + b"\x00" * 16)    # zstd frame magic
    with caplog.at_level(logging.ERROR):
        assert _run(g, tmp_path / "out") == 1
    msg = " ".join(r.getMessage() for r in caplog.records)
    assert "zstd" in msg.lower() and "decompress" in msg.lower()


def test_plain_fasta_is_not_called_compressed(tmp_path, caplog):
    g = tmp_path / "asm.fa"
    g.write_text(">c1\nACGT\n")
    with caplog.at_level(logging.ERROR):
        _run(g, tmp_path / "out")                    # may fail later (no database here); must not say "compressed"
    assert "decompress" not in " ".join(r.getMessage() for r in caplog.records).lower()
