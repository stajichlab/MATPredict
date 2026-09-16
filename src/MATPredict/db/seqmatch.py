"""Percent identity / coverage scoring between a curated accession and a re-fetched sequence."""
from __future__ import annotations

from dataclasses import dataclass

from Bio.Align import PairwiseAligner

PASS_IDENTITY_THRESHOLD = 98.0
WARN_IDENTITY_THRESHOLD = 90.0


@dataclass(frozen=True)
class MatchScore:
    """Result of aligning a query sequence against a reference sequence."""

    percent_identity: float
    coverage: float
    status: str


def _aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    return aligner


def score_match(query: str, reference: str) -> MatchScore:
    """Score a query protein sequence against a reference, tri-state pass/warn/fail."""
    if not query or not reference:
        return MatchScore(percent_identity=0.0, coverage=0.0, status="fail")

    alignment = _aligner().align(query, reference)[0]
    aligned_query, aligned_ref = str(alignment[0]), str(alignment[1])
    matches = sum(1 for a, b in zip(aligned_query, aligned_ref) if a == b and a != "-")
    aligned_columns = sum(1 for a, b in zip(aligned_query, aligned_ref) if a != "-" and b != "-")
    percent_identity = 100.0 * matches / aligned_columns if aligned_columns else 0.0
    coverage = 100.0 * aligned_columns / max(len(query), len(reference))

    if percent_identity >= PASS_IDENTITY_THRESHOLD:
        status = "pass"
    elif percent_identity >= WARN_IDENTITY_THRESHOLD:
        status = "warn"
    else:
        status = "fail"
    return MatchScore(percent_identity=round(percent_identity, 2), coverage=round(coverage, 2), status=status)
