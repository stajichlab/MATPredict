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
import shutil
import subprocess

import pytest

from MATPredict.detect import provenance
from MATPredict.detect.provenance import code_stamp_for

needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")


def _git(cwd, *args):
    subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-C", str(cwd), *args], check=True, capture_output=True)


def _project(tmp_path, version="9.9.9"):
    """A tiny checkout laid out like this project: src/MATPredict, pyproject.toml, one commit tagged v1.0."""
    pkg = tmp_path / "src" / "MATPredict"
    pkg.mkdir(parents=True)
    (pkg / "__init__.py").write_text("")
    (tmp_path / "pyproject.toml").write_text(f'[project]\nname = "MATPredict"\nversion = "{version}"\n')
    _git(tmp_path, "init", "-q")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "first")
    _git(tmp_path, "tag", "v1.0")
    return pkg


@needs_git
def test_the_stamp_names_the_tag_the_commit_and_the_source_version(tmp_path):
    pkg = _project(tmp_path, "9.9.9")
    stamp = code_stamp_for(pkg)
    assert stamp["git_describe"] == "v1.0" and stamp["dirty"] is False
    assert len(stamp["commit"]) == 40 and stamp["source_version"] == "9.9.9"


@needs_git
def test_a_commit_past_the_tag_and_a_dirty_tree_are_visible(tmp_path):
    pkg = _project(tmp_path)
    (tmp_path / "x.txt").write_text("x")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "second")
    stamp = code_stamp_for(pkg)
    assert stamp["git_describe"].startswith("v1.0-1-g") and stamp["dirty"] is False
    (tmp_path / "pyproject.toml").write_text('[project]\nversion = "0.0.1"\n')   # a tracked file changes
    stamp = code_stamp_for(pkg)
    assert stamp["git_describe"].endswith("-dirty") and stamp["dirty"] is True
    assert stamp["source_version"] == "0.0.1"          # the version of the tree that actually ran


@needs_git
def test_a_package_outside_a_git_checkout_gives_no_stamp(tmp_path):
    pkg = tmp_path / "site-packages" / "MATPredict"
    pkg.mkdir(parents=True)
    assert code_stamp_for(pkg) == {}


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


def test_run_block_carries_the_stamp_when_there_is_one_and_omits_it_otherwise(tmp_path, monkeypatch):
    db = tmp_path / "db"
    db.mkdir()
    kw = dict(clock=RunClock(), genome=None, proteins=None, db_root=db, taxonomy_source="t", sample="s",
              organism=None, taxid=None, phylum=None, parameters={})
    monkeypatch.setattr(provenance, "code_stamp", lambda: {"git_describe": "cleanrun-x", "commit": "a" * 40, "dirty": False,
                                                           "source_version": "0.7.0"})
    doc = run_block(**kw)
    assert doc["code"]["git_describe"] == "cleanrun-x" and doc["matpredict_version"] == provenance.__version__
    monkeypatch.setattr(provenance, "code_stamp", lambda: {})
    assert "code" not in run_block(**kw)


def test_the_html_provenance_shows_the_stamp():
    from MATPredict.report.html import _code_note

    assert _code_note(None) == "" and _code_note({}) == ""
    note = _code_note({"git_describe": "cleanrun-x", "commit": "0123456789abcdef"})
    assert "cleanrun-x" in note and "0123456789" in note and "abcdef" not in note
