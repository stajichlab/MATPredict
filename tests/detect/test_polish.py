from __future__ import annotations

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.polish import (
    ExonSpan, PolishModel, boundaries_agree,
    STATUS_AGREE, STATUS_DISAGREE, STATUS_SINGLE, STATUS_UNPOLISHED, classify,
)

KEY = FamilyKey("Basidiomycota", "aLocus")


def _model(exons, method="miniprot"):
    return PolishModel(
        gene_name="mfa1", family_key=KEY, role="core_MAT", contig="c1",
        start=exons[0].start, end=exons[-1].end, strand="+", exons=exons,
        identity=95.0, reference_record_id="rec1", method=method,
    )


def test_boundaries_agree_within_tolerance():
    a = _model([ExonSpan(100, 200), ExonSpan(250, 300)], method="miniprot")
    b = _model([ExonSpan(103, 198), ExonSpan(252, 297)], method="exonerate_refine")
    assert boundaries_agree(a, b, tolerance_bp=10) is True


def test_boundaries_disagree_beyond_tolerance():
    a = _model([ExonSpan(100, 200)], method="miniprot")
    b = _model([ExonSpan(150, 260)], method="exonerate_refine")
    assert boundaries_agree(a, b, tolerance_bp=10) is False


def test_boundaries_disagree_on_different_exon_count():
    a = _model([ExonSpan(100, 300)], method="miniprot")
    b = _model([ExonSpan(100, 200), ExonSpan(250, 300)], method="exonerate_refine")
    assert boundaries_agree(a, b) is False


def _model_single(start, end, method):
    """Helper for classify tests -- creates a model with a single exon."""
    return PolishModel(
        gene_name="mfa1", family_key=KEY, role="core_MAT", contig="c1",
        start=start, end=end, strand="+", exons=[ExonSpan(start, end)],
        identity=95.0, reference_record_id="rec1", method=method,
    )


def test_classify_agree_prefers_exonerate_as_canonical():
    ex = _model_single(100, 400, "exonerate_refine")
    mp = _model_single(102, 398, "miniprot")
    outcome = classify(ex, mp, tolerance_bp=10)
    assert outcome.status == STATUS_AGREE
    assert outcome.canonical is ex
    assert outcome.miniprot_model is mp


def test_classify_disagree_keeps_both_models():
    ex = _model_single(100, 400, "exonerate_refine")
    mp = _model_single(500, 900, "miniprot")
    outcome = classify(ex, mp, tolerance_bp=10)
    assert outcome.status == STATUS_DISAGREE
    assert outcome.canonical is ex
    assert outcome.exonerate_model is ex and outcome.miniprot_model is mp


def test_classify_single_tool_uses_whichever_succeeded():
    mp = _model_single(100, 400, "miniprot")
    outcome = classify(None, mp, tolerance_bp=10)
    assert outcome.status == STATUS_SINGLE
    assert outcome.canonical is mp


def test_classify_unpolished_when_neither_tool_succeeds():
    outcome = classify(None, None, tolerance_bp=10)
    assert outcome.status == STATUS_UNPOLISHED
    assert outcome.canonical is None
