import csv
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


def _acquire_via_ncbi_only(taxids, out_dir, runner, **kwargs):
    """Call `acquire_genomes` with the local-library check pointed at paths that
    are guaranteed not to exist, so these NCBI-path tests are isolated from
    whatever real BFD library/manifest may or may not be mounted on the machine
    actually running the test suite.
    """
    return acquire_genomes(
        taxids,
        out_dir,
        runner=runner,
        local_library_root=out_dir / "no_such_library",
        local_manifest_path=out_dir / "no_such_manifest.csv",
        **kwargs,
    )


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
    genomes = _acquire_via_ncbi_only([5501], tmp_path, runner)
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
    genomes = _acquire_via_ncbi_only([999999], tmp_path, runner, failures=failures)
    assert genomes == []
    assert len(failures) == 1
    assert failures[0].taxid == 999999
    assert "No assemblies found" in failures[0].reason


def test_mixed_batch_returns_successes_and_records_failures_without_raising(tmp_path):
    runner = _fake_runner(succeed_for={5501})
    failures: list[AcquisitionFailure] = []
    genomes = _acquire_via_ncbi_only([5501, 999999], tmp_path, runner, failures=failures)
    assert len(genomes) == 1
    assert genomes[0].taxid == 5501
    assert len(failures) == 1
    assert failures[0].taxid == 999999


def _write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "ASMID", "SPECIES_IN", "STRAIN", "BIOPROJECT", "NCBI_TAXONID",
        "BUSCO_LINEAGE", "PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER",
        "FAMILY", "GENUS", "SPECIES", "TRANSL_TABLE", "LOCUSTAG",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def _manifest_row(asmid: str, taxid: int, species: str) -> dict[str, str]:
    return {
        "ASMID": asmid,
        "SPECIES_IN": f"{species} test-strain",
        "NCBI_TAXONID": str(taxid),
        "SPECIES": species,
    }


def test_local_library_hit_is_used_and_ncbi_is_never_invoked(tmp_path):
    library_root = tmp_path / "library"
    library_root.mkdir()
    manifest_path = tmp_path / "manifest.csv"
    asmid = "GCA_004115165.2_Cimm211_ragoo"
    _write_manifest(manifest_path, [_manifest_row(asmid, 5501, "Coccidioides immitis")])
    fasta_path = library_root / f"{asmid}.fa.gz"
    fasta_path.write_bytes(b"fake gzip genome bytes")

    def failing_runner(cmd: list[str]) -> "subprocess.CompletedProcess[str]":
        raise AssertionError("NCBI fallback should not be invoked on a local-library hit")

    genomes = acquire_genomes(
        [5501],
        tmp_path / "ncbi_out",
        runner=failing_runner,
        local_library_root=library_root,
        local_manifest_path=manifest_path,
    )
    assert genomes == [
        AcquiredGenome(
            taxid=5501,
            accession="GCA_004115165.2",
            fasta_path=fasta_path,
            species="Coccidioides immitis",
        )
    ]


def test_local_library_prefers_unmasked_and_gcf_over_gca(tmp_path):
    library_root = tmp_path / "library"
    library_root.mkdir()
    manifest_path = tmp_path / "manifest.csv"
    gca_asmid = "GCA_000001.1_gca_asm"
    gcf_asmid = "GCF_000002.1_gcf_asm"
    _write_manifest(
        manifest_path,
        [
            _manifest_row(gca_asmid, 5501, "Coccidioides immitis"),
            _manifest_row(gcf_asmid, 5501, "Coccidioides immitis"),
        ],
    )
    # Both ASMIDs have a masked+unmasked pair; only the GCF unmasked file
    # should end up chosen (RefSeq-preference tiebreak, then unmasked
    # preference over soft-masked for the same chosen ASMID).
    (library_root / f"{gca_asmid}.fa.gz").write_bytes(b"gca unmasked")
    (library_root / f"{gcf_asmid}.masked.fasta.gz").write_bytes(b"gcf masked")
    gcf_unmasked = library_root / f"{gcf_asmid}.fa.gz"
    gcf_unmasked.write_bytes(b"gcf unmasked")

    genomes = acquire_genomes(
        [5501],
        tmp_path / "ncbi_out",
        runner=lambda cmd: (_ for _ in ()).throw(AssertionError("should not reach NCBI")),
        local_library_root=library_root,
        local_manifest_path=manifest_path,
    )
    assert len(genomes) == 1
    assert genomes[0].accession == "GCF_000002.1"
    assert genomes[0].fasta_path == gcf_unmasked


def test_local_library_miss_falls_through_to_ncbi(tmp_path):
    library_root = tmp_path / "library"
    library_root.mkdir()
    manifest_path = tmp_path / "manifest.csv"
    # Manifest exists but has no row for this taxid.
    _write_manifest(manifest_path, [_manifest_row("GCA_1.1_other", 111111, "Other species")])

    runner = _fake_runner(succeed_for={5501})
    genomes = acquire_genomes(
        [5501],
        tmp_path / "ncbi_out",
        runner=runner,
        local_library_root=library_root,
        local_manifest_path=manifest_path,
    )
    assert len(genomes) == 1
    assert genomes[0].taxid == 5501
    assert genomes[0].accession == ASSEMBLY_ACCESSION


def test_missing_manifest_path_falls_through_to_ncbi_without_raising(tmp_path):
    runner = _fake_runner(succeed_for={5501})
    genomes = acquire_genomes(
        [5501],
        tmp_path / "ncbi_out",
        runner=runner,
        local_library_root=tmp_path / "does_not_exist_library",
        local_manifest_path=tmp_path / "does_not_exist_manifest.csv",
    )
    assert len(genomes) == 1
    assert genomes[0].taxid == 5501
    assert genomes[0].accession == ASSEMBLY_ACCESSION


def test_missing_assembly_report_is_a_failure_not_a_crash(tmp_path):
    def runner(cmd: list[str]) -> "subprocess.CompletedProcess[str]":
        filename_index = cmd.index("--filename") + 1
        zip_path = Path(cmd[filename_index])
        with zipfile.ZipFile(zip_path, "w") as zf:
            zf.writestr("ncbi_dataset/data/README.md", "no report here")
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

    failures: list[AcquisitionFailure] = []
    genomes = _acquire_via_ncbi_only([5501], tmp_path, runner, failures=failures)
    assert genomes == []
    assert len(failures) == 1
    assert "assembly_data_report.jsonl" in failures[0].reason
