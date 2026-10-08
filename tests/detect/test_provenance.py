"""The `run` block of detection_report.yaml (`detect.provenance`)."""
from __future__ import annotations

import gzip

import yaml

from MATPredict.detect.pipeline import DetectionOutcome
from MATPredict.detect.provenance import RunClock, database_digest, genome_stats, run_block
from MATPredict.detect.report import write_detection_report


def test_genome_stats_counts_contigs_length_and_n50(tmp_path):
    fa = tmp_path / "g.fna"
    fa.write_text(">a\nACGTACGTAC\n>b\nACGT\nACGT\n>c\nAC\n")
    st = genome_stats(fa)
    assert (st["file"], st["contigs"], st["length_bp"], st["n50"]) == ("g.fna", 3, 20, 10)
    assert len(st["sha256"]) == 64
    gz = tmp_path / "g.fna.gz"
    gz.write_bytes(gzip.compress(fa.read_bytes()))
    assert genome_stats(gz)["length_bp"] == 20


def test_database_digest_ignores_candidates_and_tracks_content(tmp_path):
    db = tmp_path / "db"
    (db / "Asco" / "r1").mkdir(parents=True)
    (db / "candidates").mkdir()
    (db / "Asco" / "r1" / "metadata.yaml").write_text("a: 1\n")
    first = database_digest(db)
    (db / "candidates" / "x.yaml").write_text("new candidate\n")
    assert database_digest(db) == first
    assert first["records"] == 1
    (db / "Asco" / "r1" / "metadata.yaml").write_text("a: 2\n")
    assert database_digest(db)["content_sha256"] != first["content_sha256"]


def test_run_block_survives_an_unreadable_genome(tmp_path):
    db = tmp_path / "db"
    db.mkdir()
    doc = run_block(clock=RunClock(), genome=tmp_path / "missing.fna", proteins=None, db_root=db,
                    taxonomy_source="test", sample=None, organism=None, taxid=4837, phylum="Mucoromycota",
                    parameters={})
    assert doc["sample"] == "missing"
    assert "error" in doc["genome"]


def test_run_is_the_first_key_and_optional(tmp_path):
    out = tmp_path / "r.yaml"
    write_detection_report(DetectionOutcome(results=[]), out, run={"sample": "s1"})
    doc = yaml.safe_load(out.read_text())
    assert next(iter(doc)) == "run" and doc["run"]["sample"] == "s1"
    write_detection_report(DetectionOutcome(results=[]), out)
    assert "run" not in yaml.safe_load(out.read_text())
