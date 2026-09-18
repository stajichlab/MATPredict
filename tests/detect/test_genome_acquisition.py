import json
import subprocess
import zipfile
from pathlib import Path

from MATPredict.detect.genome_acquisition import (
    AcquiredGenome,
    AcquisitionFailure,
    acquire_genomes,
)

ASSEMBLY_ACCESSION = "GCA_000149405.2"
SPECIES_NAME = "Coccidioides immitis"
GENOMIC_FASTA_TEXT = ">CM000000.1 test contig\nACGTACGTACGT\n"


def _write_fake_datasets_zip(zip_path: Path, accession: str, species: str) -> None:
    """Build a zip with the same internal layout the real `datasets` CLI produces."""
    report = {"accession": accession, "organism": {"organismName": species}}
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("ncbi_dataset/data/assembly_data_report.jsonl", json.dumps(report) + "\n")
        zf.writestr(
            f"ncbi_dataset/data/{accession}/{accession}_ASM_genomic.fna",
            GENOMIC_FASTA_TEXT,
        )


def _fake_runner(succeed_for: set[int] | None = None, accession: str = ASSEMBLY_ACCESSION,
                  species: str = SPECIES_NAME):
    """Fake subprocess runner: mimics `datasets download`'s side effect of writing a
    zip file to the path given by `--filename`, instead of making a real network call
    or spawning a real process -- following this repo's convention (see
    `tests/db/test_ncbi_client.py`'s `_fake_transport`) of faking at the transport
    boundary rather than mocking library internals.
    """

    def runner(cmd: list[str]) -> "subprocess.CompletedProcess[str]":
        taxid = int(cmd[cmd.index("taxon") + 1])
        filename_index = cmd.index("--filename") + 1
        zip_path = Path(cmd[filename_index])
        if succeed_for is not None and taxid not in succeed_for:
            return subprocess.CompletedProcess(cmd, returncode=1, stdout="", stderr="No assemblies found")
        _write_fake_datasets_zip(zip_path, accession, species)
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

    return runner


def test_successful_acquisition_returns_right_acquired_genome(tmp_path):
    runner = _fake_runner()
    genomes = acquire_genomes([5501], tmp_path, runner=runner)
    assert genomes == [
        AcquiredGenome(
            taxid=5501,
            accession=ASSEMBLY_ACCESSION,
            fasta_path=tmp_path / "5501" / "extracted" / "ncbi_dataset" / "data"
            / ASSEMBLY_ACCESSION / f"{ASSEMBLY_ACCESSION}_ASM_genomic.fna",
            species=SPECIES_NAME,
        )
    ]
    assert genomes[0].fasta_path.read_text() == GENOMIC_FASTA_TEXT


def test_unresolvable_taxid_is_skipped_and_reported_as_failure(tmp_path):
    runner = _fake_runner(succeed_for=set())
    failures: list[AcquisitionFailure] = []
    genomes = acquire_genomes([999999], tmp_path, runner=runner, failures=failures)
    assert genomes == []
    assert len(failures) == 1
    assert failures[0].taxid == 999999
    assert "No assemblies found" in failures[0].reason


def test_mixed_batch_returns_successes_and_records_failures_without_raising(tmp_path):
    runner = _fake_runner(succeed_for={5501})
    failures: list[AcquisitionFailure] = []
    genomes = acquire_genomes([5501, 999999], tmp_path, runner=runner, failures=failures)
    assert len(genomes) == 1
    assert genomes[0].taxid == 5501
    assert len(failures) == 1
    assert failures[0].taxid == 999999


def test_missing_assembly_report_is_a_failure_not_a_crash(tmp_path):
    def runner(cmd: list[str]) -> "subprocess.CompletedProcess[str]":
        filename_index = cmd.index("--filename") + 1
        zip_path = Path(cmd[filename_index])
        with zipfile.ZipFile(zip_path, "w") as zf:
            zf.writestr("ncbi_dataset/data/README.md", "no report here")
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

    failures: list[AcquisitionFailure] = []
    genomes = acquire_genomes([5501], tmp_path, runner=runner, failures=failures)
    assert genomes == []
    assert len(failures) == 1
    assert "assembly_data_report.jsonl" in failures[0].reason
