# tests/db/test_seqmatch_coverage.py
"""A perfect identity over a sliver of the protein is not a pass.

`score_match` decided pass/warn/fail on percent_identity ALONE and ignored the
coverage it had just computed. That let a badly wrong record validate clean:
during the 2026-09-20 ingest of OR965930.1, exons were recorded in genomic
order instead of transcript order, which corrupts every minus-strand
multi-exon gene. The result validated as

    algA  identity=100.0  coverage=1.30   pass
    glrA  identity=100.0  coverage=0.86   pass

-- 100% identical over less than one percent of the protein, because a single
exon happened to land in frame. Had that been accepted, proteins.faa would
have carried garbage translations for the two genes being added.

It is a WARNING, not a failure. Low coverage is legitimate in real curated
records: Aspergillus nidulans 162425_fgsc-a4_MAT_combined has genes at 3.0%,
18.9% and 30.5% coverage with 100% identity, where the curated CDS span
deliberately covers part of a larger deposited protein. Failing those would
invalidate correct curation; warning on them surfaces the question without
overruling a curator.
"""
from __future__ import annotations

from MATPredict.db.seqmatch import WARN_COVERAGE_THRESHOLD, score_match


def test_a_full_length_perfect_match_passes():
    seq = "MKVLAAGIVGSTQWDFPYRNCLEHM" * 4
    result = score_match(query=seq, reference=seq)
    assert result.percent_identity == 100.0
    assert result.coverage == 100.0
    assert result.status == "pass"


def test_perfect_identity_over_a_sliver_is_warned_not_passed():
    # The real shape of the exon-order bug: a short in-frame fragment matching
    # perfectly against a long deposited protein.
    reference = "MKVLAAGIVGSTQWDFPYRNCLEHM" * 20   # 500 aa
    query = reference[:12]                           # 12 aa, 2.4% coverage
    result = score_match(query=query, reference=reference)
    assert result.percent_identity == 100.0
    assert result.coverage < WARN_COVERAGE_THRESHOLD
    assert result.status == "warn", "identity alone must not carry a pass"


def test_a_low_identity_match_still_fails_regardless_of_coverage():
    # Coverage must not rescue a genuinely wrong sequence either.
    reference = "MKVLAAGIVGSTQWDFPYRNCLEHM" * 4
    query = "PQRSTVWYACDEFGHIKLMNPQRST" * 4
    result = score_match(query=query, reference=reference)
    assert result.status == "fail"


def test_coverage_at_the_threshold_still_passes():
    reference = "MKVLAAGIVGSTQWDFPYRNCLEHM" * 20
    keep = int(len(reference) * WARN_COVERAGE_THRESHOLD / 100.0) + 1
    result = score_match(query=reference[:keep], reference=reference)
    assert result.coverage >= WARN_COVERAGE_THRESHOLD
    assert result.status == "pass"
