from __future__ import annotations

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.polish import ExonSpan, PolishModel, boundaries_agree

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
