"""The HTML report shows where a locus's own genes lie and how far the span runs beyond them (field `core_span`)."""
from __future__ import annotations

from pathlib import Path

import yaml

from MATPredict.report.html import render_genome_report
from tests.report.test_report_html import _parse

FIXTURES = Path(__file__).parent / "fixtures"


def _text(core_span) -> str:
    doc = yaml.safe_load((FIXTURES / "basidio_hd_pr.yaml").read_text())
    doc["detected"][0]["core_span"] = core_span
    return " ".join(_parse(render_genome_report(doc, sample="x")).text)


def test_a_span_beyond_the_genes_is_shown_in_plain_words():
    # T48-F: genes 118,043-130,301; the span runs 22,895 bp past them
    text = _text({"start": 118043, "end": 130301, "beyond_core_bp": 22895})
    assert "Own genes" in text and "118,043" in text and "130,301" in text
    assert "22.9 kb" in text and "beyond" in text


def test_a_span_equal_to_its_genes_shows_the_genes_without_a_warning():
    text = _text({"start": 100, "end": 200, "beyond_core_bp": 0})
    assert "Own genes" in text and "beyond" not in text


def test_reports_from_before_the_field_do_not_show_it():
    assert "Own genes" not in _text(None)


def _text_sup(supported) -> str:
    doc = yaml.safe_load((FIXTURES / "basidio_hd_pr.yaml").read_text())
    doc["detected"][0]["supported_span"] = supported
    return " ".join(_parse(render_genome_report(doc, sample="x")).text)


def test_a_span_beyond_its_supported_extent_is_shown_with_the_floor():
    # T48-F-like: 22,895 bp of the span carried only by hits below 33 bits
    text = _text_sup({"start": 118043, "end": 130301, "min_bitscore": 33.0, "beyond_supported_bp": 22895})
    assert "Supported span" in text and "118,043" in text and "130,301" in text
    assert "22.9 kb" in text and "33 bits" in text


def test_a_fully_supported_span_is_not_called_out():
    text = _text_sup({"start": 100, "end": 200, "min_bitscore": 33.0, "beyond_supported_bp": 0})
    assert "Supported span" not in text


def test_reports_without_the_field_do_not_show_a_supported_span():
    assert "Supported span" not in _text_sup(None)
