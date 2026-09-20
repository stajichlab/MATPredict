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

    def fake_run_pipeline(genome_fasta, proteome_fasta, taxid, db_root, reference_fasta, evidence_diagnostics_path=None):
        # The decompressed fasta must exist and be a real, uncompressed FASTA
        # file inside $SCRATCH at call time.
        assert genome_fasta.exists()
        assert str(genome_fasta).startswith(str(scratch))
        assert genome_fasta.read_text().startswith(">contig1")
        # Task 4: run_batch wires evidence_diagnostics_path to
        # genome_out_dir/evidence_diagnostics.jsonl for every genome.
        assert evidence_diagnostics_path is not None
        assert evidence_diagnostics_path.name == "evidence_diagnostics.jsonl"
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

    def flaky_run_pipeline(genome_fasta, proteome_fasta, taxid, db_root, reference_fasta, evidence_diagnostics_path=None):
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


def test_run_batch_continues_after_one_genome_fails_to_decompress(tmp_path, monkeypatch):
    """A genuinely corrupt/truncated .gz source must not abort the rest of the
    batch, and must be recorded via the same `failures` mechanism a
    run_pipeline failure uses -- not a second, separate reporting path.
    """
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setenv("SCRATCH", str(scratch))

    corrupt = tmp_path / "corrupt.fa.gz"
    corrupt.write_bytes(b"this is not gzip data at all")
    good = tmp_path / "good.fa.gz"
    _write_gz_fasta(good)
    genomes = [
        _genome(1, "CORRUPT", corrupt),
        _genome(2, "GOOD", good),
    ]

    calls = []

    def fake_run_pipeline(genome_fasta, proteome_fasta, taxid, db_root, reference_fasta, evidence_diagnostics_path=None):
        calls.append(taxid)
        return DetectionOutcome(results=[])

    out_dir = tmp_path / "out"
    failures: list[GenomeRunFailure] = []
    run_batch(
        genomes,
        tmp_path / "db",
        tmp_path / "reference.faa",
        out_dir,
        run_pipeline=fake_run_pipeline,
        failures=failures,
    )

    # run_pipeline was never called for the corrupt genome, but processing
    # continued and the good genome still ran and produced a report.
    assert calls == [2]
    assert (out_dir / "2_GOOD" / "detected_loci.gff3").exists()

    # The decompression failure was recorded the same way a run_pipeline
    # failure would be -- one GenomeRunFailure in the same `failures` list.
    assert len(failures) == 1
    assert failures[0].taxid == 1
    assert failures[0].accession == "CORRUPT"

    # No partial/leftover decompressed file left behind on scratch.
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


def _run_batch_recording_gff3(tmp_path, monkeypatch, **run_batch_kwargs):
    """Run a one-genome batch with `write_detection_gff3` replaced by a
    recorder, and return the keyword arguments it received."""
    import MATPredict.detect.batch_runner as batch_runner

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setenv("SCRATCH", str(scratch))

    genome_path = tmp_path / "genome.fa.gz"
    _write_gz_fasta(genome_path)
    genome = _genome(7, "GCA_000000007.1", genome_path)

    recorded: list[dict] = []

    def fake_write_gff3(outcome, out_path, **kwargs):
        recorded.append(kwargs)

    monkeypatch.setattr(batch_runner, "write_detection_gff3", fake_write_gff3)

    run_batch(
        [genome], tmp_path / "db", tmp_path / "reference.faa", tmp_path / "out",
        run_pipeline=lambda **kwargs: DetectionOutcome(results=[]),
        **run_batch_kwargs,
    )
    return recorded, scratch


def test_run_batch_omits_genome_fasta_by_default(tmp_path, monkeypatch):
    """`emit_cds_fasta` defaults False so every existing rollout caller keeps
    producing byte-identical output and pays no extra genome parse."""
    recorded, _ = _run_batch_recording_gff3(tmp_path, monkeypatch)
    assert recorded == [{}]


def test_run_batch_passes_decompressed_genome_when_emit_cds_fasta(tmp_path, monkeypatch):
    """When asked for CDS output, `run_batch` must reuse the copy it already
    decompressed into `$SCRATCH` for `run_pipeline` -- re-decompressing or
    re-fetching the genome to get a second readable copy would double the
    I/O for no benefit, and the `.gz` source is not what `run_pipeline`
    itself was handed."""
    recorded, scratch = _run_batch_recording_gff3(
        tmp_path, monkeypatch, emit_cds_fasta=True,
    )
    assert list(recorded[0]) == ["genome_fasta"]
    passed = recorded[0]["genome_fasta"]
    assert str(passed).startswith(str(scratch))
    assert not str(passed).endswith(".gz")


def test_run_batch_emit_cds_fasta_is_keyword_only(tmp_path, monkeypatch):
    """Positional passing must not be possible: `run_batch`'s existing
    positional arguments are load-bearing at several call sites, and a new
    positional flag could silently bind to the wrong one."""
    import inspect
    parameter = inspect.signature(run_batch).parameters["emit_cds_fasta"]
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default is False


# ---------------------------------------------------------------------------
# Task 4: the batch path must be able to pay the same narrowed cost as the
# single-genome CLI path (`matpredict detect --phylum`), and must be able to
# choose its own evidence floor.
# ---------------------------------------------------------------------------


def _mini_db(root: Path) -> Path:
    """A two-phylum curated database just large enough for routing and for
    `build_reference_fasta`'s family filter to have something to filter.

    Written per-test rather than pointing at the live `db/` because these
    tests assert WHICH proteins reach the query set, and the live database is
    under active curation -- an assertion on its contents would break on the
    next accepted record. The live-database measurement lives in
    `test_reference_fasta.py`, stated as a relationship for the same reason.
    """
    for phylum, order, record_id, locus, scope, protein in (
        ("Mucoromycota", "Mucorales", "muco_rec", "sexMP", 4827, "MMM"),
        ("Ascomycota", "Saccharomycetales", "asco_rec", "MATsc", 4932, "AAA"),
    ):
        (root / phylum).mkdir(parents=True, exist_ok=True)
        (root / phylum / "order.yml").write_text(
            f"phylum: {phylum}\n"
            "loci:\n"
            f"  - locus_name: {locus}\n"
            "    vocabulary_type: idiomorph\n"
            "    idiomorph_values: [a, b]\n"
            "    genes: [g1]\n"
            f"    taxonomic_scope: [{scope}]\n"
        )
        record_dir = root / phylum / order / record_id
        record_dir.mkdir(parents=True, exist_ok=True)
        (record_dir / "metadata.yaml").write_text(
            f"record_id: {record_id}\nmating_type:\n  locus_name: {locus}\n"
        )
        (record_dir / "proteins.faa").write_text(
            f">{record_id}|gene_index=0|name=g1|role=core_MAT\n{protein}\n"
        )
    return root


def _two_genomes(tmp_path: Path) -> list[AcquiredGenome]:
    genomes = []
    for i, taxid in enumerate((4827, 4932), start=1):
        path = tmp_path / f"batch_genome{i}.fa.gz"
        _write_gz_fasta(path)
        genomes.append(_genome(taxid, f"GCA_00000000{i}.1", path))
    return genomes


def test_run_batch_phylum_routes_every_genome_to_that_phylum(tmp_path, monkeypatch):
    """Requirement 1: `--phylum` on a batch must restrict ROUTING for every
    genome in it, not just the first, and must not depend on each genome's
    taxid -- the second genome's taxid here is an Ascomycota one, and it must
    still be routed to Mucoromycota because the operator said so.
    """
    monkeypatch.setenv("SCRATCH", str(tmp_path / "scratch"))
    (tmp_path / "scratch").mkdir()
    db_root = _mini_db(tmp_path / "db")

    seen = []

    def fake_run_pipeline(*, routing, **kwargs):
        seen.append(routing)
        return DetectionOutcome(results=[])

    run_batch(
        _two_genomes(tmp_path), db_root, tmp_path / "out" / "_reference.faa",
        tmp_path / "out", run_pipeline=fake_run_pipeline, phylum="Mucoromycota",
    )

    assert len(seen) == 2
    for routing in seen:
        assert routing.routing_mode == "explicit_phylum"
        assert {f.key.phylum for f in routing.families} == {"Mucoromycota"}
    # Routed ONCE and reused, not re-derived per genome.
    assert seen[0] is seen[1]


def test_run_batch_phylum_builds_the_restricted_reference_fasta_once(tmp_path, monkeypatch):
    """Requirement 1 + the binding build-once constraint: the query set must
    hold only the requested phylum's curated proteins, and must be built one
    time for the whole batch -- rebuilding it per genome would walk every
    `db/*/*/*/proteins.faa` 813 times in a full rollout.
    """
    import MATPredict.detect.batch_runner as batch_runner

    monkeypatch.setenv("SCRATCH", str(tmp_path / "scratch"))
    (tmp_path / "scratch").mkdir()
    db_root = _mini_db(tmp_path / "db")

    real_build = batch_runner.build_reference_fasta
    builds = []

    def counting_build(*args, **kwargs):
        builds.append((args, kwargs))
        return real_build(*args, **kwargs)

    monkeypatch.setattr(batch_runner, "build_reference_fasta", counting_build)

    references = []

    def fake_run_pipeline(*, reference_fasta, **kwargs):
        references.append(reference_fasta)
        return DetectionOutcome(results=[])

    reference_path = tmp_path / "out" / "_reference.faa"
    run_batch(
        _two_genomes(tmp_path), db_root, reference_path, tmp_path / "out",
        run_pipeline=fake_run_pipeline, phylum="Mucoromycota",
    )

    assert len(builds) == 1, "the reference FASTA must be built once per batch"
    # The phylum goes into the destination filename, so a restricted query set
    # can never overwrite the name an unrestricted one uses, and two array
    # tasks scoped to different phyla cannot clobber each other's query set.
    restricted_path = tmp_path / "out" / "_reference.Mucoromycota.faa"
    assert references == [restricted_path, restricted_path]
    assert not reference_path.exists()
    text = restricted_path.read_text()
    assert ">muco_rec|gene0|g1" in text
    assert "asco_rec" not in text


def test_run_batch_without_phylum_self_routes_unchanged(tmp_path, monkeypatch):
    """Requirement 4: a mixed-phylum batch must still work. With no `phylum`,
    `run_batch` passes NO routing and NO reference FASTA of its own -- it hands
    `run_pipeline` exactly what the caller built, and `run_pipeline` self-routes
    per genome as before. The fake below accepts no `routing` keyword on
    purpose: passing one would raise TypeError and fail this test.
    """
    import MATPredict.detect.batch_runner as batch_runner

    monkeypatch.setenv("SCRATCH", str(tmp_path / "scratch"))
    (tmp_path / "scratch").mkdir()
    db_root = _mini_db(tmp_path / "db")

    monkeypatch.setattr(
        batch_runner, "build_reference_fasta",
        lambda *a, **k: pytest.fail("run_batch must not rebuild the reference FASTA"),
    )

    references = []

    def fake_run_pipeline(genome_fasta, proteome_fasta, taxid, db_root, reference_fasta,
                          evidence_diagnostics_path=None):
        references.append(reference_fasta)
        return DetectionOutcome(results=[])

    caller_reference = tmp_path / "caller_reference.faa"
    caller_reference.write_text(">rec|gene0|g1\nMMM\n")
    run_batch(
        _two_genomes(tmp_path), db_root, caller_reference, tmp_path / "out",
        run_pipeline=fake_run_pipeline,
    )

    assert references == [caller_reference, caller_reference]


def test_run_batch_rejects_unknown_phylum_before_any_genome_runs(tmp_path, monkeypatch):
    """A `phylum` matching no curated family is a usage error, exactly as in
    `cli._cmd_detect`. It must be raised BEFORE the per-genome loop: inside it,
    the loop's own try/except (which exists so one bad genome cannot sink the
    batch) would swallow it once per genome and an 813-genome job would run to
    completion having searched nothing.
    """
    monkeypatch.setenv("SCRATCH", str(tmp_path / "scratch"))
    (tmp_path / "scratch").mkdir()
    db_root = _mini_db(tmp_path / "db")

    failures: list[GenomeRunFailure] = []

    def fake_run_pipeline(**kwargs):
        pytest.fail("no genome may be processed under an unroutable phylum")

    with pytest.raises(ValueError, match="Chytridiomycota"):
        run_batch(
            _two_genomes(tmp_path), db_root, tmp_path / "out" / "_reference.faa",
            tmp_path / "out", run_pipeline=fake_run_pipeline, failures=failures,
            phylum="Chytridiomycota",
        )

    assert failures == []


def test_run_batch_threads_evidence_floor_to_run_pipeline(tmp_path, monkeypatch):
    """Requirement 2: a batch must be able to tighten or loosen the floor."""
    from MATPredict.detect.pipeline import EvidenceFloor

    monkeypatch.setenv("SCRATCH", str(tmp_path / "scratch"))
    (tmp_path / "scratch").mkdir()

    floor = EvidenceFloor(min_hits=1, require_core_role=False)
    seen = []

    def fake_run_pipeline(*, evidence_floor, **kwargs):
        seen.append(evidence_floor)
        return DetectionOutcome(results=[])

    run_batch(
        _two_genomes(tmp_path), tmp_path / "db", tmp_path / "reference.faa",
        tmp_path / "out", run_pipeline=fake_run_pipeline, evidence_floor=floor,
    )

    assert seen == [floor, floor]


def test_run_batch_omits_evidence_floor_when_none(tmp_path, monkeypatch):
    """`evidence_floor=None` must mean "let `run_pipeline` use its own
    default", which is enforced by NOT passing the keyword at all. Building an
    `EvidenceFloor()` here instead would create a second copy of the default
    that can silently drift from `run_pipeline`'s -- the exact bug Task 3 hit
    with argparse defaults. The fake below accepts no `evidence_floor`
    keyword, so passing one fails this test.
    """
    monkeypatch.setenv("SCRATCH", str(tmp_path / "scratch"))
    (tmp_path / "scratch").mkdir()

    calls = []

    def fake_run_pipeline(genome_fasta, proteome_fasta, taxid, db_root, reference_fasta,
                          evidence_diagnostics_path=None):
        calls.append(taxid)
        return DetectionOutcome(results=[])

    run_batch(
        _two_genomes(tmp_path), tmp_path / "db", tmp_path / "reference.faa",
        tmp_path / "out", run_pipeline=fake_run_pipeline,
    )

    assert len(calls) == 2


def test_run_batch_new_scope_parameters_are_keyword_only():
    """Same reasoning as `emit_cds_fasta`: `run_batch`'s positionals are
    load-bearing, and a new positional could silently bind to the wrong one.
    Both default to today's behaviour.
    """
    import inspect

    parameters = inspect.signature(run_batch).parameters
    for name in ("phylum", "evidence_floor"):
        assert parameters[name].kind is inspect.Parameter.KEYWORD_ONLY
        assert parameters[name].default is None


def test_run_batch_phylum_still_isolates_one_bad_genome(tmp_path, monkeypatch):
    """The per-genome try/except isolation is binding and must survive the new
    phylum path: a genome that fails under `--phylum` is recorded and the batch
    continues.
    """
    monkeypatch.setenv("SCRATCH", str(tmp_path / "scratch"))
    (tmp_path / "scratch").mkdir()
    db_root = _mini_db(tmp_path / "db")

    genomes = _two_genomes(tmp_path)
    ok = []

    def fake_run_pipeline(*, taxid, **kwargs):
        if taxid == genomes[0].taxid:
            raise RuntimeError("boom")
        ok.append(taxid)
        return DetectionOutcome(results=[])

    failures: list[GenomeRunFailure] = []
    run_batch(
        genomes, db_root, tmp_path / "out" / "_reference.faa", tmp_path / "out",
        run_pipeline=fake_run_pipeline, failures=failures, phylum="Mucoromycota",
    )

    assert ok == [genomes[1].taxid]
    assert [f.taxid for f in failures] == [genomes[0].taxid]
