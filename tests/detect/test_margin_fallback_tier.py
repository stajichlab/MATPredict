"""A thin vote margin lowers a `high` call to `medium`; an identity gap used as a stand-in margin does not (#137)."""
from MATPredict.detect.pipeline import _thin_vote_margin_caps_tier as caps


def test_a_thin_vote_margin_lowers_high():
    assert caps(2.0, False, 5.0, "high")


def test_a_wide_vote_margin_keeps_high():
    assert not caps(5.0, False, 5.0, "high")


def test_the_resolution_identity_gap_never_lowers_the_tier():
    assert not caps(0.8, True, 5.0, "high")


def test_only_high_is_lowered_and_no_margin_means_no_cap():
    assert not caps(2.0, False, 5.0, "medium")
    assert not caps(None, False, 5.0, "high")
