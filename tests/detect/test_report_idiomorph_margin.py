# tests/detect/test_report_idiomorph_margin.py
"""The report carries how close the idiomorph call was, and what it resolved.

The curator's ruling was to make the call, cap the tier when it is close, and
ALSO report the ambiguity -- explicitly because hybrids and novel loci are
things worth discovering, and a resolution that silently discarded the losing
evidence would hide exactly the cases worth a second look. A number only in a
diagnostics file nobody opens does not satisfy that.
"""
from __future__ import annotations

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.idiomorph import IdiomorphResolution
from MATPredict.detect.pipeline import DetectionResult
from MATPredict.detect.report import _result_doc

KEY = FamilyKey("Mucoromycota", "MAT")


def _result(**kw) -> DetectionResult:
    base = dict(
        family_key=KEY, contig="c1", start=1, end=100, confidence="high",
        idiomorph="Plus", ambiguous_with=[], genes_found=["tptA", "sexP"],
        genes_missing=[], fragmented=False,
    )
    base.update(kw)
    return DetectionResult(**base)


def test_a_locus_needing_no_resolution_reports_a_null_margin():
    doc = _result_doc(_result())
    assert doc["idiomorph_margin"] is None
    assert doc["idiomorph_resolutions"] == []


def test_the_margin_and_both_members_reach_the_report():
    resolution = IdiomorphResolution(
        contig="c1", winner="sexM", loser="sexP",
        winner_identity=30.60, loser_identity=28.26,
        overlap_fraction=1.0, winner_coverage=None, loser_coverage=23.6,
    )
    doc = _result_doc(_result(
        idiomorph="Minus", idiomorph_margin=resolution.margin,
        idiomorph_resolutions=[resolution],
    ))
    assert abs(doc["idiomorph_margin"] - 2.34) < 1e-9
    assert doc["idiomorph_resolutions"] == [{
        "contig": "c1",
        "winner": "sexM",
        "loser": "sexP",
        "winner_identity": 30.60,
        "loser_identity": 28.26,
        "overlap_fraction": 1.0,
        "winner_coverage": None,
        "loser_coverage": 23.6,
    }]
