"""Task 4: the SLURM entry point must expose the batch path's new scope.

`scripts/run_detection_batch.py` is what an 813-genome rollout actually runs,
so `run_batch`'s `phylum` and `evidence_floor` parameters are only reachable in
practice if this driver exposes them. The flag names are asserted to be
IDENTICAL to `matpredict detect`'s: an operator comparing a single-genome
debug run against a batch run must not have to translate flag spellings, and a
silently different spelling here would be discovered only mid-rollout.

The script is loaded by path because `scripts/` is not an importable package.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script():
    spec = importlib.util.spec_from_file_location(
        "run_detection_batch_script", _REPO_ROOT / "scripts" / "run_detection_batch.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _genomes_json(tmp_path: Path) -> Path:
    path = tmp_path / "batch.json"
    path.write_text(
        json.dumps(
            [
                {
                    "taxid": 4827,
                    "accession": "GCA_000000001.1",
                    "fasta_path": str(tmp_path / "g1.fa.gz"),
                    "species": "sp1",
                }
            ]
        )
    )
    return path


def _run(module, tmp_path, monkeypatch, extra_argv: list[str]):
    """Run the script's `main` with `run_batch` and `build_reference_fasta`
    replaced by recorders, and return (run_batch kwargs, build call count)."""
    recorded: dict = {}
    builds: list[tuple] = []

    def fake_run_batch(genomes, db_root, reference_fasta, out_dir, **kwargs):
        recorded["genomes"] = genomes
        recorded["db_root"] = db_root
        recorded["reference_fasta"] = reference_fasta
        recorded["out_dir"] = out_dir
        recorded.update(kwargs)

    def fake_build(db_root, out_path, **kwargs):
        builds.append((db_root, out_path, kwargs))
        return out_path

    monkeypatch.setattr(module, "run_batch", fake_run_batch)
    monkeypatch.setattr(module, "build_reference_fasta", fake_build)

    argv = [
        "--genomes-json", str(_genomes_json(tmp_path)),
        "--db-root", str(tmp_path / "db"),
        "--out-dir", str(tmp_path / "out"),
        *extra_argv,
    ]
    module.main(argv)
    return recorded, builds


def test_script_flag_names_match_matpredict_detect():
    """Requirement 3: same spellings as the single-genome CLI."""
    module = _load_script()
    parser = module._build_parser()
    options = {
        option for action in parser._actions for option in action.option_strings
    }
    for flag in ("--phylum", "--min-hits", "--min-identity",
                 "--require-core-role", "--no-require-core-role"):
        assert flag in options, f"{flag} is missing from the batch driver"


def test_script_passes_phylum_through_and_does_not_prebuild_reference(tmp_path, monkeypatch):
    """With `--phylum`, `run_batch` builds the ONE restricted reference FASTA
    itself (it is the component that knows the routed families). The driver
    must therefore NOT also build an unrestricted one first -- that would walk
    every curated `proteins.faa` for a file nothing ever reads.
    """
    module = _load_script()
    recorded, builds = _run(module, tmp_path, monkeypatch, ["--phylum", "Mucoromycota"])

    assert recorded["phylum"] == "Mucoromycota"
    assert builds == []
    assert recorded["reference_fasta"] == tmp_path / "out" / "_reference.faa"


def test_script_without_phylum_builds_the_unrestricted_reference_once(tmp_path, monkeypatch):
    """Requirement 4: a mixed-phylum batch keeps working exactly as before --
    one unrestricted reference FASTA, built once by the driver, and no phylum
    forced on `run_batch`.
    """
    module = _load_script()
    recorded, builds = _run(module, tmp_path, monkeypatch, [])

    assert recorded["phylum"] is None
    assert len(builds) == 1
    assert builds[0][2] == {}


def test_script_evidence_floor_defaults_come_from_evidence_floor(tmp_path, monkeypatch):
    """The driver must not restate the floor's defaults. If it did, the batch
    rollout would keep running the OLD permissive floor the moment
    `EvidenceFloor`'s defaults changed -- which is exactly how Task 3's
    defaults reached this path unnoticed in the first place.
    """
    from MATPredict.detect.pipeline import EvidenceFloor

    module = _load_script()
    recorded, _ = _run(module, tmp_path, monkeypatch, [])

    assert recorded["evidence_floor"] == EvidenceFloor()


def test_script_evidence_floor_flags_override_the_defaults(tmp_path, monkeypatch):
    from MATPredict.detect.pipeline import EvidenceFloor

    module = _load_script()
    recorded, _ = _run(
        module, tmp_path, monkeypatch,
        ["--min-hits", "1", "--min-identity", "30", "--no-require-core-role"],
    )

    assert recorded["evidence_floor"] == EvidenceFloor(
        min_hits=1, min_identity=30.0, require_core_role=False
    )


def test_script_require_core_role_can_be_switched_back_on(tmp_path, monkeypatch):
    module = _load_script()
    recorded, _ = _run(module, tmp_path, monkeypatch, ["--require-core-role"])
    assert recorded["evidence_floor"].require_core_role is True
