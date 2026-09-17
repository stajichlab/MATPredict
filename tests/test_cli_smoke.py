from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import duckdb
import yaml

from MATPredict.__main__ import main
from MATPredict.db.validate import validate_record


def test_detect_subcommand_registered():
    from MATPredict.__main__ import build_parser
    parser = build_parser()
    args = parser.parse_args(["detect", "--genome", "g.fa", "--out-dir", "/tmp/x"])
    assert args.command == "detect"
    assert args.genome == "g.fa"


def test_help_exits_zero(capsys):
    exit_code = main(["curate-db", "--help"])
    assert exit_code == 0
    captured = capsys.readouterr()
    # Verify action-specific content from curate-db subparser (not just top-level help)
    # Actions are registered in db/cli.py: propose, validate, accept, reject, build-gff, build-duckdb, release
    assert "propose" in captured.out or "validate" in captured.out


def test_build_duckdb_subcommand_runs_against_real_db(tmp_path, monkeypatch):
    monkeypatch.setenv("MATPREDICT_DB_ROOT", str(Path(__file__).resolve().parents[1] / "db"))
    out_path = tmp_path / "matpredict.duckdb"
    exit_code = main(["curate-db", "build-duckdb", "--out", str(out_path)])
    assert exit_code == 0
    assert out_path.exists()
    con = duckdb.connect(str(out_path))
    con.execute("SELECT 1 FROM locus_record LIMIT 1")  # doesn't raise, table exists
    con.close()


# --- Ruling 1: validate must force validation.status to needs_review on check failure ---

_ESUMMARY_SUPPRESSED = (
    '{"result": {"uids": ["1"], "1": {"status": "suppressed"}}}'
)


def _fake_requests_get(monkeypatch, response_text: str):
    class _FakeResponse:
        text = response_text
        status_code = 200

    def _get(url):
        return _FakeResponse()

    monkeypatch.setattr("MATPredict.db.cli.requests.get", _get)


def _fake_taxonkit_run(cmd, **kwargs):
    """Stand-in for subprocess.run(["taxonkit", ...]) — never touch the real binary/taxdump."""
    return SimpleNamespace(returncode=0, stdout="4837\tk__Fungi;p__Mucoromycota;...;s__Phycomyces_blakesleeanus\n")


def _patch_taxonomy_runner(monkeypatch):
    """Replace validate_record's bound default `taxonomy_runner=subprocess.run` with a fake.

    `_cmd_validate` calls `validate_record(record, ncbi=ncbi, uniprot=uniprot)` with no
    `taxonomy_runner` argument, so it always falls through to validate.py's default value,
    which was bound to the real `subprocess.run` function object at module-import time.
    Monkeypatching `subprocess.run` (globally or via `MATPredict.db.taxonomy.subprocess.run`)
    has *no effect* here: a function's default argument value is evaluated once, at
    definition time, and is not re-looked-up on each call — patching the `subprocess`
    module's `run` attribute afterward does not change what's already stored in
    `validate_record.__defaults__`. Patching `__defaults__` itself is the only way to
    intercept this without touching validate.py's public interface; pytest's monkeypatch
    fixture restores the original tuple automatically at teardown.
    """
    monkeypatch.setattr(validate_record, "__defaults__", (_fake_taxonkit_run,))


def _write_candidate(tmp_path, phylum: str, record_id: str, validation_status: str) -> Path:
    record = {
        "record_id": record_id,
        "taxonomy": {"taxid": 4837},
        "locus": {
            "coordinate_provenance": "published_explicit",
            "core": {
                "segments": [
                    {"segment_index": 0, "sequence_source": {"type": "insdc_nucleotide", "accession": "EU009461.1"}}
                ],
            },
        },
        "genes": [],
        "validation": {"status": validation_status, "rejection_reason": None},
    }
    candidate_dir = tmp_path / "db" / "candidates" / phylum / record_id
    candidate_dir.mkdir(parents=True)
    (candidate_dir / "metadata.yaml").write_text(yaml.safe_dump(record, sort_keys=False))
    return candidate_dir / "metadata.yaml"


def test_validate_forces_needs_review_on_accession_not_resolved(tmp_path, monkeypatch):
    monkeypatch.setenv("MATPREDICT_DB_ROOT", str(tmp_path / "db"))
    monkeypatch.setenv("MATPREDICT_CACHE_DIR", str(tmp_path / "cache"))
    metadata_path = _write_candidate(tmp_path, "Mucoromycota", "rec1", validation_status="accepted")
    _fake_requests_get(monkeypatch, _ESUMMARY_SUPPRESSED)
    _patch_taxonomy_runner(monkeypatch)

    exit_code = main(["curate-db", "validate", "--phylum", "Mucoromycota", "--record-id", "rec1"])
    assert exit_code == 0

    record = yaml.safe_load(metadata_path.read_text())
    assert record["validation"]["accession_resolved"] is False
    assert record["validation"]["status"] == "needs_review"


def test_validate_never_un_rejects_a_rejected_record(tmp_path, monkeypatch):
    monkeypatch.setenv("MATPREDICT_DB_ROOT", str(tmp_path / "db"))
    monkeypatch.setenv("MATPREDICT_CACHE_DIR", str(tmp_path / "cache"))
    metadata_path = _write_candidate(tmp_path, "Mucoromycota", "rec2", validation_status="rejected")
    _fake_requests_get(monkeypatch, _ESUMMARY_SUPPRESSED)
    _patch_taxonomy_runner(monkeypatch)

    exit_code = main(["curate-db", "validate", "--phylum", "Mucoromycota", "--record-id", "rec2"])
    assert exit_code == 0

    record = yaml.safe_load(metadata_path.read_text())
    assert record["validation"]["accession_resolved"] is False
    assert record["validation"]["status"] == "rejected"


# --- Ruling 2: build-gff wires GFF3/GenBank/FASTA export for an accepted record ---

_EFETCH_FASTA = ">AAB12345.1\nMKTAYIAKQRQISFVKSHFSRQ\n"


def test_build_gff_subcommand_writes_all_three_files(tmp_path, monkeypatch):
    monkeypatch.setenv("MATPREDICT_DB_ROOT", str(tmp_path / "db"))
    monkeypatch.setenv("MATPREDICT_CACHE_DIR", str(tmp_path / "cache"))

    record = {
        "record_id": "rec3",
        "locus": {
            "core": {
                "segments": [
                    {"segment_index": 0, "sequence_source": {"seq_region": "scaffold_3"}, "start": 100, "end": 500}
                ],
            },
        },
        "genes": [
            {
                "gene_index": 0,
                "name": "sexP",
                "role": "core_MAT",
                "present": True,
                "segment_index": 0,
                "start": 150,
                "end": 250,
                "strand": "+",
                "protein_accession": "ncbi_protein:AAB12345.1",
            },
        ],
        "validation": {"status": "accepted", "rejection_reason": None},
    }
    record_dir = tmp_path / "db" / "Mucoromycota" / "Mucorales" / "rec3"
    record_dir.mkdir(parents=True)
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(record, sort_keys=False))

    class _FakeResponse:
        text = _EFETCH_FASTA
        status_code = 200

    monkeypatch.setattr("MATPredict.db.cli.requests.get", lambda url: _FakeResponse())

    exit_code = main([
        "curate-db", "build-gff",
        "--phylum", "Mucoromycota",
        "--order-or-family", "Mucorales",
        "--record-id", "rec3",
    ])
    assert exit_code == 0

    gff3_path = record_dir / "locus.gff3"
    gbk_path = record_dir / "locus.gbk"
    proteins_path = record_dir / "proteins.faa"
    assert gff3_path.exists()
    assert gbk_path.exists()
    assert proteins_path.exists()

    assert "sexP" in gff3_path.read_text()
    assert "sexP" in gbk_path.read_text()
    assert "MKTAYIAKQRQISFVKSHFSRQ" in proteins_path.read_text()
