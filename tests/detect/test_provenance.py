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


# --- the code stamp (issue #139): the reported version is the installed metadata, so the checkout is read directly ---
import os
import shutil
import subprocess

import pytest

from MATPredict.detect import provenance
from MATPredict.detect.provenance import code_stamp_for

needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")


@pytest.fixture(autouse=True)
def _hermetic_git(monkeypatch, tmp_path):
    """The host's git configuration (signing, hooks, safe.directory, global ignores) must not reach these tests, and git
    must not look for a repository above the temporary directory."""
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", os.devnull)
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path.parent))


def _git(cwd, *args):
    subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false",
                    "-c", "tag.gpgsign=false", "-c", "core.hooksPath=/dev/null", "-C", str(cwd), *args],
                   check=True, capture_output=True)


def _project(tmp_path, version="9.9.9", tag="v1.0"):
    """A tiny checkout laid out like this project: src/MATPredict, pyproject.toml, one commit (tagged unless tag=None)."""
    pkg = tmp_path / "src" / "MATPredict"
    pkg.mkdir(parents=True)
    (pkg / "__init__.py").write_text("")
    (tmp_path / "pyproject.toml").write_text(f'[project]\nname = "MATPredict"\nversion = "{version}"\n')
    (tmp_path / "db").mkdir()
    (tmp_path / "db" / "r.yml").write_text("a: 1\n")
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs" / "n.md").write_text("notes\n")
    _git(tmp_path, "init", "-q")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "first")
    if tag:
        _git(tmp_path, "tag", tag)
    return pkg


@needs_git
def test_the_stamp_names_the_tag_the_commit_and_the_source_version(tmp_path):
    pkg = _project(tmp_path, "9.9.9")
    stamp = code_stamp_for(pkg)
    assert stamp["git_describe"] == "v1.0" and stamp["dirty"] is False
    assert len(stamp["commit"]) == 40 and stamp["source_version"] == "9.9.9"
    head = subprocess.run(["git", "-C", str(tmp_path), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    assert stamp["commit"] == head


@needs_git
def test_a_commit_past_the_tag_names_the_nearest_tag_and_the_distance(tmp_path):
    pkg = _project(tmp_path)
    (tmp_path / "x.txt").write_text("x")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "second")
    stamp = code_stamp_for(pkg)
    assert stamp["git_describe"].startswith("v1.0-1-g") and stamp["dirty"] is False


@needs_git
def test_dirty_means_code_changed_including_untracked_files_and_not_a_touched_doc_or_result(tmp_path):
    pkg = _project(tmp_path)
    assert code_stamp_for(pkg)["dirty"] is False
    (tmp_path / "docs" / "n.md").write_text("changed notes\n")                 # a tracked file outside the code paths
    (tmp_path / "results").mkdir()
    (tmp_path / "results" / "summary.tsv").write_text("untracked result\n")
    assert code_stamp_for(pkg)["dirty"] is False
    (pkg / "dropped_in.py").write_text("x = 1\n")                              # an untracked module changes behaviour
    assert code_stamp_for(pkg)["dirty"] is True
    (pkg / "dropped_in.py").unlink()
    assert code_stamp_for(pkg)["dirty"] is False
    (tmp_path / "db" / "r.yml").write_text("a: 2\n")                           # an edited database file is code too
    assert code_stamp_for(pkg)["dirty"] is True
    _git(tmp_path, "checkout", "--", "db/r.yml")
    (tmp_path / "pyproject.toml").write_text('[project]\nname = "MATPredict"\nversion = "0.0.1"\n')
    stamp = code_stamp_for(pkg)
    assert stamp["dirty"] is True and stamp["source_version"] == "0.0.1"        # the version of the tree that actually ran


@needs_git
def test_a_repository_with_no_tags_falls_back_to_the_abbreviated_hash(tmp_path):
    pkg = _project(tmp_path, tag=None)
    stamp = code_stamp_for(pkg)
    assert stamp["commit"].startswith(stamp["git_describe"]) and len(stamp["git_describe"]) < 40


@needs_git
def test_a_detached_head_and_a_git_worktree_are_stamped(tmp_path):
    main = tmp_path / "main"
    main.mkdir()
    _project(main)
    wt = tmp_path / "wt"
    _git(main, "worktree", "add", "-q", "--detach", str(wt), "v1.0")            # `.git` is a file here, HEAD is detached
    assert (wt / ".git").is_file()
    stamp = code_stamp_for(wt / "src" / "MATPredict")
    assert stamp["git_describe"] == "v1.0" and stamp["dirty"] is False and stamp["source_version"] == "9.9.9"


@needs_git
def test_a_package_outside_any_git_checkout_gives_no_stamp(tmp_path):
    pkg = tmp_path / "site-packages" / "MATPredict"
    pkg.mkdir(parents=True)
    assert code_stamp_for(pkg) == {}                  # GIT_CEILING_DIRECTORIES keeps git from finding an outer repository


@needs_git
def test_a_repository_around_somebody_elses_venv_is_not_this_projects_checkout(tmp_path):
    _git(tmp_path, "init", "-q")                       # a repo, but the package is not at src/MATPredict
    pkg = tmp_path / ".venv" / "lib" / "site-packages" / "MATPredict"
    pkg.mkdir(parents=True)
    assert code_stamp_for(pkg) == {}


def test_a_missing_git_binary_gives_no_stamp_and_no_error(tmp_path, monkeypatch):
    def boom(*a, **k):
        raise FileNotFoundError("git")

    monkeypatch.setattr(provenance.subprocess, "run", boom)
    assert code_stamp_for(tmp_path) == {}


def test_an_unreadable_checkout_is_logged_not_silent(tmp_path, monkeypatch, caplog):
    root = tmp_path / "proj"
    (root / "src" / "MATPredict" / "detect").mkdir(parents=True)
    (root / ".git").mkdir()                                            # looks like a checkout; git cannot read it
    fake = root / "src" / "MATPredict" / "detect" / "provenance.py"
    fake.write_text("")
    monkeypatch.setattr(provenance, "__file__", str(fake))
    monkeypatch.setattr(provenance, "code_stamp_for", lambda package_dir: {})
    provenance.code_stamp.cache_clear()
    try:
        with caplog.at_level("WARNING"):
            assert provenance.code_stamp() == {}
        assert any("could not read the git identity" in r.getMessage() for r in caplog.records)
    finally:
        provenance.code_stamp.cache_clear()


@needs_git
def test_the_stamp_of_the_development_checkout_matches_git(monkeypatch):
    provenance.code_stamp.cache_clear()
    stamp = provenance.code_stamp()
    provenance.code_stamp.cache_clear()
    if not stamp:
        pytest.skip("the package is not running from a git checkout")
    top = subprocess.run(["git", "-C", str(provenance.Path(provenance.__file__).parent), "rev-parse", "--show-toplevel"],
                         capture_output=True, text=True).stdout.strip()
    head = subprocess.run(["git", "-C", top, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    assert stamp["commit"] == head and stamp["source_version"]


def test_run_block_carries_the_stamp_when_there_is_one_and_omits_it_otherwise(tmp_path, monkeypatch):
    db = tmp_path / "db"
    db.mkdir()
    kw = dict(clock=RunClock(), genome=None, proteins=None, db_root=db, taxonomy_source="t", sample="s",
              organism=None, taxid=None, phylum=None, parameters={})
    monkeypatch.setattr(provenance, "code_stamp", lambda: {"git_describe": "cleanrun-x", "commit": "a" * 40, "dirty": False,
                                                           "source_version": "0.7.0"})
    doc = run_block(**kw)
    assert doc["code"]["git_describe"] == "cleanrun-x"
    monkeypatch.setattr(provenance, "code_stamp", lambda: {})
    assert "code" not in run_block(**kw)


def test_run_block_records_the_running_version(tmp_path, monkeypatch):
    db = tmp_path / "db"
    db.mkdir()
    kw = dict(clock=RunClock(), genome=None, proteins=None, db_root=db, taxonomy_source="t", sample="s",
              organism=None, taxid=None, phylum=None, parameters={})
    monkeypatch.setattr(provenance, "running_version", lambda: "7.7.7")
    assert run_block(**kw)["matpredict_version"] == "7.7.7"


def _pyproject(root, name="MATPredict", version="0.7.0"):
    pkg = root / "src" / "MATPredict"
    pkg.mkdir(parents=True)
    (root / "pyproject.toml").write_text(f'[project]\nname = "{name}"\nversion = "{version}"\n')
    return pkg


def test_the_source_version_is_the_checkouts_pyproject_version_and_needs_no_git(tmp_path):
    assert provenance.source_version_for(_pyproject(tmp_path / "a")) == ("0.7.0" if provenance.tomllib else None)
    if provenance.tomllib is None:
        pytest.skip("tomllib needs Python 3.11")
    assert provenance.source_version_for(_pyproject(tmp_path / "b", version="1.2.3")) == "1.2.3"


def test_the_source_version_is_none_for_a_wheel_another_project_or_an_unreadable_file(tmp_path):
    if provenance.tomllib is None:
        pytest.skip("tomllib needs Python 3.11")
    wheel = tmp_path / "venv" / "lib" / "site-packages" / "MATPredict"
    wheel.mkdir(parents=True)
    assert provenance.source_version_for(wheel) is None                              # not a src/ layout
    assert provenance.source_version_for(_pyproject(tmp_path / "other", name="somethingelse")) is None
    pkg = _pyproject(tmp_path / "bad")
    (tmp_path / "bad" / "pyproject.toml").write_text("not [valid toml")
    assert provenance.source_version_for(pkg) is None
    (tmp_path / "bad" / "pyproject.toml").unlink()
    assert provenance.source_version_for(pkg) is None                                # no pyproject at all


def test_the_running_version_falls_back_to_the_installed_metadata(monkeypatch):
    monkeypatch.setattr(provenance, "source_version_for", lambda package_dir: None)
    provenance.running_version.cache_clear()
    try:
        assert provenance.running_version() == provenance.__version__
        monkeypatch.setattr(provenance, "source_version_for", lambda package_dir: "8.8.8")
        provenance.running_version.cache_clear()
        assert provenance.running_version() == "8.8.8"
    finally:
        provenance.running_version.cache_clear()


def test_the_html_provenance_shows_the_stamp_the_source_version_mismatch_and_a_dirty_tree():
    from MATPredict.report.html import _code_note

    assert _code_note(None) == "" and _code_note({}) == ""
    note = _code_note({"git_describe": "cleanrun-x", "commit": "0123456789abcdef"})
    assert "cleanrun-x" in note and "0123456789" in note and "abcdef" not in note
    # the #139 symptom: the recorded version is stale compared with the checkout that ran
    stale = _code_note({"git_describe": "cleanrun-x", "commit": "0123456789abcdef", "source_version": "0.6.3"}, "0.6.1")
    assert "source 0.6.3" in stale
    assert "source" not in _code_note({"git_describe": "v1", "commit": "0123456789", "source_version": "0.7.0"}, "0.7.0")
    assert "uncommitted changes" in _code_note({"git_describe": "v1", "commit": "0123456789", "dirty": True})
    # an untagged checkout describes as the abbreviated hash: shown once, not twice
    once = _code_note({"git_describe": "24829be", "commit": "24829be1a2c3d4e5f60718293a4b5c6d7e8f9012"})
    assert once.count("24829be") == 1
