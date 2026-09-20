"""Task 4: the SLURM entry point must expose the batch path's new scope.

`scripts/run_detection_batch.py` is what an 813-genome rollout actually runs,
so `run_batch`'s `phylum` and `evidence_floor` parameters are only reachable in
practice if this driver exposes them. The flag names are asserted to be
IDENTICAL to `matpredict detect`'s: an operator comparing a single-genome
debug run against a batch run must not have to translate flag spellings, and a
silently different spelling here would be discovered only mid-rollout.

`scripts/run_detection_batch.slurm` is covered too: it is the DOCUMENTED
rollout entry point, so a flag the driver accepts but the wrapper cannot
forward is, in practice, a flag the rollout does not have. Those tests run the
real wrapper with a stub `pixi` on PATH rather than reading its text, so they
assert the command it actually invokes.

The scripts are loaded/run by path because `scripts/` is not an importable
package.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

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


def _option_strings(parser) -> set[str]:
    return {option for action in parser._actions for option in action.option_strings}


def _detect_parser():
    """The REAL `matpredict detect` subparser, built the way `__main__` builds
    it. Comparing against a hardcoded list of flag names instead would compare
    the driver against a copy of itself: renaming `--min-hits` in `cli.py`
    would leave such a test green while the two entry points silently
    diverged, which is the exact drift this assertion exists to catch.
    """
    import argparse

    from MATPredict.detect.cli import register_subcommands

    parser = argparse.ArgumentParser()
    register_subcommands(parser.add_subparsers(dest="command"))
    return parser._subparsers._group_actions[0].choices["detect"]


def test_script_flag_names_match_matpredict_detect():
    """Requirement 3: the scope and evidence-floor flags must be spelled
    EXACTLY as `matpredict detect` spells them, asserted against that parser
    itself."""
    module = _load_script()
    driver_options = _option_strings(module._build_parser())
    detect_options = _option_strings(_detect_parser())

    shared = {"--phylum", "--min-hits", "--min-identity",
              "--require-core-role", "--no-require-core-role"}
    missing_from_cli = shared - detect_options
    assert not missing_from_cli, (
        f"{sorted(missing_from_cli)} is not a `matpredict detect` flag -- either "
        "the CLI renamed it (and the driver must follow) or this test names a "
        "flag that never existed"
    )
    missing_from_driver = shared - driver_options
    assert not missing_from_driver, (
        f"{sorted(missing_from_driver)} is missing from the batch driver"
    )


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


# ---------------------------------------------------------------------------
# The SLURM wrapper is the DOCUMENTED rollout entry point. If it cannot pass
# the new flags through, an operator following its own usage line runs the
# 813-genome rollout unnarrowed -- all 181 curated proteins per genome -- with
# nothing in the output saying so. These tests run the real script with a stub
# `pixi` on PATH, so they assert what the wrapper actually invokes.
# ---------------------------------------------------------------------------

_SLURM_WRAPPER = _REPO_ROOT / "scripts" / "run_detection_batch.slurm"


def _run_wrapper(tmp_path: Path, args: list[str]):
    import subprocess

    stub_bin = tmp_path / "bin"
    stub_bin.mkdir(exist_ok=True)
    recorded = tmp_path / "pixi_argv.txt"
    stub = stub_bin / "pixi"
    stub.write_text(
        "#!/usr/bin/bash\n"
        f'printf "%s\\n" "$@" > {recorded}\n'
    )
    stub.chmod(0o755)

    env = {
        "PATH": f"{stub_bin}:/usr/bin:/bin",
        "SCRATCH": str(tmp_path / "scratch"),
        "HOME": str(tmp_path),
    }
    completed = subprocess.run(
        ["/usr/bin/bash", str(_SLURM_WRAPPER), *args],
        cwd=_REPO_ROOT, env=env, capture_output=True, text=True,
    )
    argv = recorded.read_text().splitlines() if recorded.exists() else []
    return completed, argv


def test_slurm_wrapper_forwards_scope_and_floor_flags(tmp_path):
    completed, argv = _run_wrapper(
        tmp_path,
        ["batch.json", "db", "out",
         "--phylum", "Mucoromycota", "--min-hits", "1", "--no-require-core-role"],
    )

    assert completed.returncode == 0, completed.stderr
    assert argv[:3] == ["run", "python", "scripts/run_detection_batch.py"]
    assert argv[3:] == [
        "--genomes-json", "batch.json",
        "--db-root", "db",
        "--out-dir", "out",
        "--phylum", "Mucoromycota",
        "--min-hits", "1",
        "--no-require-core-role",
    ]


def test_slurm_wrapper_keeps_the_three_positional_form_working(tmp_path):
    """Existing invocations must not break."""
    completed, argv = _run_wrapper(tmp_path, ["batch.json", "db", "out"])

    assert completed.returncode == 0, completed.stderr
    assert argv[3:] == [
        "--genomes-json", "batch.json", "--db-root", "db", "--out-dir", "out",
    ]


def test_slurm_wrapper_logs_the_scope_it_ran_with(tmp_path):
    """Whether the query set was narrowed changes what the run's negatives
    mean, and the sbatch command line is not recorded in the job log. The log
    must therefore state the scope itself."""
    scoped, _ = _run_wrapper(
        tmp_path, ["batch.json", "db", "out", "--phylum", "Mucoromycota"]
    )
    assert "phylum=Mucoromycota" in scoped.stdout

    unscoped, _ = _run_wrapper(tmp_path, ["batch.json", "db", "out"])
    assert "phylum=<none>" in unscoped.stdout


def test_slurm_wrapper_logs_the_scope_for_the_equals_spelling(tmp_path):
    scoped, argv = _run_wrapper(
        tmp_path, ["batch.json", "db", "out", "--phylum=Mucoromycota"]
    )
    assert "phylum=Mucoromycota" in scoped.stdout
    assert argv[-1] == "--phylum=Mucoromycota"


def test_slurm_wrapper_still_rejects_missing_positionals(tmp_path):
    completed, argv = _run_wrapper(tmp_path, ["batch.json", "db"])
    assert completed.returncode == 1
    assert "usage:" in completed.stderr
    assert argv == []
