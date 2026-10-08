"""`detect` writes report.html by default; `--no-html` / `$MATPREDICT_HTML=0` turn it off; a report error never
fails the run."""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from MATPredict.detect.cli import _write_html_report, html_report_enabled

FIXTURE = Path(__file__).parent / "fixtures" / "phycomyces_classifier.yaml"


@pytest.mark.parametrize("flag, env, expected", [
    (None, None, True),
    (None, "0", False),
    (None, "off", False),
    (None, "1", True),
    (True, "0", True),     # the flag wins over the variable
    (False, None, False),
])
def test_html_default_flag_and_env(monkeypatch, flag, env, expected):
    if env is None:
        monkeypatch.delenv("MATPREDICT_HTML", raising=False)
    else:
        monkeypatch.setenv("MATPREDICT_HTML", env)
    assert html_report_enabled(SimpleNamespace(html=flag)) is expected


def test_report_written_next_to_the_yaml(tmp_path, monkeypatch):
    monkeypatch.delenv("MATPREDICT_HTML", raising=False)
    (tmp_path / "detection_report.yaml").write_text(FIXTURE.read_text())
    _write_html_report(SimpleNamespace(html=None, pdf=False), tmp_path)
    assert "Phycomyces blakesleeanus" in (tmp_path / "report.html").read_text()


def test_no_html_writes_nothing(tmp_path, monkeypatch):
    monkeypatch.setenv("MATPREDICT_HTML", "0")
    (tmp_path / "detection_report.yaml").write_text(FIXTURE.read_text())
    _write_html_report(SimpleNamespace(html=None, pdf=False), tmp_path)
    assert not (tmp_path / "report.html").exists()


def test_report_failure_does_not_raise(tmp_path, monkeypatch, caplog):
    monkeypatch.delenv("MATPREDICT_HTML", raising=False)
    (tmp_path / "detection_report.yaml").write_text("not: a report\n")
    _write_html_report(SimpleNamespace(html=None, pdf=False), tmp_path)
    assert not (tmp_path / "report.html").exists()
    assert "report not written" in caplog.text


def test_missing_pdf_engine_keeps_the_html(tmp_path, monkeypatch, caplog):
    from MATPredict.report import pdf
    monkeypatch.setattr(pdf, "available_engine", lambda engine=None: (_ for _ in ()).throw(pdf.NoPdfEngine("none")))
    (tmp_path / "detection_report.yaml").write_text(FIXTURE.read_text())
    _write_html_report(SimpleNamespace(html=True, pdf=True), tmp_path)
    assert (tmp_path / "report.html").exists()
    assert not (tmp_path / "report.pdf").exists()
    assert "PDF not written" in caplog.text
