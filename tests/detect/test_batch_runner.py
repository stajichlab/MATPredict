from __future__ import annotations

import gzip
import os
from pathlib import Path

import pytest

from MATPredict.detect.batch_runner import (
    GenomeRunFailure,
    estimate_seconds_from_gzipped_size,
    plan_batches,
    run_batch,
)
from MATPredict.detect.genome_acquisition import AcquiredGenome
from MATPredict.detect.pipeline import DetectionOutcome


def _fake_estimator(seconds_by_path: dict[Path, float]):
    def estimate(fasta_path: Path) -> float:
        return seconds_by_path[fasta_path]

    return estimate


def _genome(taxid: int, accession: str, fasta_path: Path, species: str = "sp") -> AcquiredGenome:
    return AcquiredGenome(taxid=taxid, accession=accession, fasta_path=fasta_path, species=species)


def test_plan_batches_groups_toward_target_runtime(tmp_path):
    a = tmp_path / "a.fasta"
    b = tmp_path / "b.fasta"
    a.write_text("x")
    b.write_text("x")
    genomes = [_genome(1, "A", a), _genome(2, "B", b)]
    estimate = _fake_estimator({a: 2700.0, b: 2700.0})

    batches = plan_batches(genomes, target_seconds_per_job=5400, estimate_seconds=estimate)

    assert len(batches) == 1
    assert len(batches[0]) == 2


def test_plan_batches_splits_when_over_target(tmp_path):
    a = tmp_path / "a.fasta"
    b = tmp_path / "b.fasta"
    a.write_text("x")
    b.write_text("x")
    genomes = [_genome(1, "A", a), _genome(2, "B", b)]
    # Each genome alone already exceeds the target -- must not be merged.
    estimate = _fake_estimator({a: 6000.0, b: 6000.0})

    batches = plan_batches(genomes, target_seconds_per_job=5400, estimate_seconds=estimate)

    assert len(batches) == 2
    assert batches[0] == [genomes[0]]
    assert batches[1] == [genomes[1]]


def test_plan_batches_preserves_order_across_batches(tmp_path):
    paths = [tmp_path / f"g{i}.fasta" for i in range(4)]
    for p in paths:
        p.write_text("x")
    genomes = [_genome(i, f"ACC{i}", p) for i, p in enumerate(paths)]
    # 2500s each -> two per batch (5000s) fits under a 5400s target.
    estimate = _fake_estimator({p: 2500.0 for p in paths})

    batches = plan_batches(genomes, target_seconds_per_job=5400, estimate_seconds=estimate)

    assert [g.taxid for batch in batches for g in batch] == [0, 1, 2, 3]
    assert [len(b) for b in batches] == [2, 2]


def test_estimate_seconds_from_gzipped_size_scales_linearly(tmp_path):
    # Reference point: ~12 MB gzipped -> 96.4s. A file at half that size should
    # estimate to roughly half the time.
    small = tmp_path / "small.fa.gz"
    small.write_bytes(b"0" * (6 * 1_000_000))
    big = tmp_path / "big.fa.gz"
    big.write_bytes(b"0" * (12 * 1_000_000))

    small_estimate = estimate_seconds_from_gzipped_size(small)
    big_estimate = estimate_seconds_from_gzipped_size(big)

    assert big_estimate == pytest.approx(96.4, rel=0.05)
    assert small_estimate == pytest.approx(big_estimate / 2, rel=0.05)


def _write_gz_fasta(path: Path, text: str = ">contig1\nACGT\n") -> None:
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def test_run_batch_decompresses_calls_pipeline_and_writes_reports(tmp_path, monkeypatch):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setenv("SCRATCH", str(scratch))

    genome_a = tmp_path / "genomeA.fa.gz"
    genome_b = tmp_path / "genomeB.masked.fasta.gz"
    _write_gz_fasta(genome_a)
    _write_gz_fasta(genome_b)
    genomes = [
        _genome(111, "GCA_000000001.1", genome_a, species="sp1"),
        _genome(222, "GCA_000000002.2", genome_b, species="sp2"),
    ]

    calls = []

    def fake_run_pipeline(genome_fasta, proteome_fasta, taxid, db_root, reference_fasta):
        # The decompressed fasta must exist and be a real, uncompressed FASTA
        # file inside $SCRATCH at call time.
        assert genome_fasta.exists()
        assert str(genome_fasta).startswith(str(scratch))
        assert genome_fasta.read_text().startswith(">contig1")
        calls.append((taxid, genome_fasta))
        return DetectionOutcome(results=[])

    out_dir = tmp_path / "out"
    db_root = tmp_path / "db"
    reference_fasta = tmp_path / "reference.faa"
    reference_fasta.write_text("")

    run_batch(genomes, db_root, reference_fasta, out_dir, run_pipeline=fake_run_pipeline)

    assert len(calls) == 2
    assert {t for t, _ in calls} == {111, 222}

    for genome in genomes:
        genome_dir = out_dir / f"{genome.taxid}_{genome.accession}"
        assert (genome_dir / "detected_loci.gff3").exists()
        assert (genome_dir / "detection_report.yaml").exists()

    # Decompressed copies must be cleaned up after each genome's run.
    leftover = list(scratch.iterdir())
    assert leftover == []


def test_run_batch_continues_after_one_genome_fails(tmp_path, monkeypatch):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setenv("SCRATCH", str(scratch))

    genome_a = tmp_path / "bad.fa.gz"
    genome_b = tmp_path / "good.fa.gz"
    _write_gz_fasta(genome_a)
    _write_gz_fasta(genome_b)
    genomes = [
        _genome(1, "BAD", genome_a),
        _genome(2, "GOOD", genome_b),
    ]

    def flaky_run_pipeline(genome_fasta, proteome_fasta, taxid, db_root, reference_fasta):
        if taxid == 1:
            raise RuntimeError("boom")
        return DetectionOutcome(results=[])

    out_dir = tmp_path / "out"
    failures: list[GenomeRunFailure] = []
    run_batch(
        genomes,
        tmp_path / "db",
        tmp_path / "reference.faa",
        out_dir,
        run_pipeline=flaky_run_pipeline,
        failures=failures,
    )

    assert len(failures) == 1
    assert failures[0].taxid == 1
    assert failures[0].accession == "BAD"
    # The second genome's report was still written.
    assert (out_dir / "2_GOOD" / "detected_loci.gff3").exists()
    # Scratch is still cleaned up even for the failing genome.
    assert list(scratch.iterdir()) == []


def test_run_batch_raises_loudly_when_scratch_unset(tmp_path, monkeypatch):
    monkeypatch.delenv("SCRATCH", raising=False)
    genome = _genome(1, "A", tmp_path / "a.fa.gz")

    with pytest.raises(RuntimeError, match="SCRATCH"):
        run_batch([genome], tmp_path / "db", tmp_path / "reference.faa", tmp_path / "out")


def test_run_batch_sanitizes_accession_with_dots(tmp_path, monkeypatch):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setenv("SCRATCH", str(scratch))

    genome_path = tmp_path / "genome.fa.gz"
    _write_gz_fasta(genome_path)
    genome = _genome(5334, "GCA_004115165.2", genome_path)

    def fake_run_pipeline(**kwargs):
        return DetectionOutcome(results=[])

    out_dir = tmp_path / "out"
    run_batch([genome], tmp_path / "db", tmp_path / "reference.faa", out_dir, run_pipeline=fake_run_pipeline)

    assert (out_dir / "5334_GCA_004115165.2" / "detected_loci.gff3").exists()
