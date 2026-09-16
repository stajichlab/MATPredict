from __future__ import annotations

from MATPredict.db.seqmatch import score_match


def test_identical_sequences_pass():
    seq = "MKTAYIAKQRQISFVKSHFSRQ"
    result = score_match(query=seq, reference=seq)
    assert result.status == "pass"
    assert result.percent_identity == 100.0
    assert result.coverage == 100.0


def test_single_mismatch_still_passes_at_high_identity():
    reference = "MKTAYIAKQRQISFVKSHFSRQ"
    query = "MKTAYIAKQRQISFVKSHFSRA"  # 1 of 22 differ => ~95.5% identity
    result = score_match(query=query, reference=reference)
    assert result.status in {"pass", "warn"}
    assert result.percent_identity > 90.0


def test_unrelated_sequences_fail():
    result = score_match(query="AAAAAAAAAAAAAAAAAAAAAA", reference="MKTAYIAKQRQISFVKSHFSRQ")
    assert result.status == "fail"


def test_empty_query_or_reference_fails():
    assert score_match(query="", reference="MKTAYIAKQRQ").status == "fail"
    assert score_match(query="MKTAYIAKQRQ", reference="").status == "fail"
