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


def test_detect_evidence_floor_flags_parse_with_defaults():
    from MATPredict.__main__ import build_parser
    parser = build_parser()
    args = parser.parse_args(["detect", "--genome", "g.fa", "--out-dir", "/tmp/x"])
    assert args.evidence_diagnostics is None
    assert args.min_hits == 1
    assert args.min_identity is None
    assert args.require_core_role is False


def test_detect_evidence_floor_flags_parse_when_given():
    from MATPredict.__main__ import build_parser
    parser = build_parser()
    args = parser.parse_args([
        "detect", "--genome", "g.fa", "--out-dir", "/tmp/x",
        "--evidence-diagnostics", "/tmp/x/diag.jsonl",
        "--min-hits", "3", "--min-identity", "45.5", "--require-core-role",
    ])
    assert args.evidence_diagnostics == "/tmp/x/diag.jsonl"
    assert args.min_hits == 3
    assert args.min_identity == 45.5
    assert args.require_core_role is True


def test_detect_default_evidence_floor_args_reconstruct_the_no_op_default():
    """When none of the new flags are passed, `_cmd_detect`'s
    `EvidenceFloor(min_hits=args.min_hits, min_identity=args.min_identity,
    require_core_role=args.require_core_role)` must equal `EvidenceFloor()`
    exactly -- i.e. a caller that doesn't pass the new flags gets identical
    `run_pipeline` behavior to before this change."""
    from MATPredict.__main__ import build_parser
    from MATPredict.detect.pipeline import EvidenceFloor
    parser = build_parser()
    args = parser.parse_args(["detect", "--genome", "g.fa", "--out-dir", "/tmp/x"])
    reconstructed = EvidenceFloor(
        min_hits=args.min_hits, min_identity=args.min_identity,
        require_core_role=args.require_core_role,
    )
    assert reconstructed == EvidenceFloor()
    assert args.evidence_diagnostics is None


def test_detect_missing_required_args_exits_code_1():
    """Regression test: --genome/--out-dir moved from argparse required=True to
    app-level check so detect benchmark could share the same parser. This changed
    the exit code for missing args from argparse's standard 2 to 1 (generic
    exception handler). Lock in this behavior so future refactors can't silently
    flip it again without test coverage."""
    exit_code = main(["detect"])
    assert exit_code == 1


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


def test_validate_forces_needs_review_on_undeclared_gene_name(tmp_path, monkeypatch, capsys):
    """A curated gene whose name is not in its locus's order.yml `genes` list is
    silently dropped by detect.search._attribute, so it can never be found in any
    genome. `validate` must refuse to let it through review unnoticed."""
    monkeypatch.setenv("MATPREDICT_DB_ROOT", str(tmp_path / "db"))
    monkeypatch.setenv("MATPREDICT_CACHE_DIR", str(tmp_path / "cache"))
    metadata_path = _write_candidate(tmp_path, "Mucoromycota", "rec9", validation_status="accepted")
    record = yaml.safe_load(metadata_path.read_text())
    record["locus"]["coordinate_provenance"] = "not_available"  # skip the network-dependent checks
    record["mating_type"] = {"locus_name": "MAT", "idiomorphs": ["Plus"], "system": "heterothallic"}
    record["genes"] = [{"gene_index": 0, "name": "sexP_undeclared", "role": "core_MAT", "present": True}]
    metadata_path.write_text(yaml.safe_dump(record, sort_keys=False))

    order_path = tmp_path / "db" / "Mucoromycota" / "order.yml"
    order_path.parent.mkdir(parents=True)
    order_path.write_text(yaml.safe_dump({
        "phylum": "Mucoromycota",
        "loci": [{
            "locus_name": "MAT",
            "vocabulary_type": "enum",
            "idiomorph_values": ["Plus", "Minus"],
            "taxonomic_scope": [4827],
            "genes": [{"name": "sexP", "role": "core_MAT"}],
        }],
    }))
    _patch_taxonomy_runner(monkeypatch)

    assert main(["curate-db", "validate", "--phylum", "Mucoromycota", "--record-id", "rec9"]) == 0
    assert "sexP_undeclared" in capsys.readouterr().out
    assert yaml.safe_load(metadata_path.read_text())["validation"]["status"] == "needs_review"


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
