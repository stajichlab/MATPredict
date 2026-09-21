"""Percent identity / coverage scoring between a curated accession and a re-fetched sequence."""
from __future__ import annotations

from dataclasses import dataclass

from Bio.Align import PairwiseAligner

PASS_IDENTITY_THRESHOLD = 98.0
WARN_IDENTITY_THRESHOLD = 90.0

WARN_COVERAGE_THRESHOLD = 80.0
"""Below this coverage a perfect identity is downgraded to a warning.

80.0 rather than something tighter because low coverage is legitimate in real
curated records: across the 171 accepted genes carrying a coverage value, 164
sit at 95% or above, but the remainder run 3.0-86.7% and are correct -- the
curated CDS span deliberately covers part of a larger deposited protein
(Aspergillus nidulans 162425_fgsc-a4_MAT_combined, Schizophyllum commune
5270_521_bLocus_b1). 80.0 leaves the highest of those (86.68%) passing while
still catching the failure this exists for: a 0.86% coverage match that
validated as a clean pass.
"""


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

    # Identity alone cannot carry a pass. A perfect match over a sliver of the
    # protein means the record's coordinates do not describe the protein it
    # claims: during the 2026-09-20 ingest of OR965930.1, exons recorded in
    # genomic instead of transcript order produced
    # `identity=100.0 coverage=0.86 pass` for glrA, because one exon happened
    # to land in frame. Accepting that would have written garbage translations
    # into proteins.faa for the very genes being added.
    #
    # A WARNING, never a failure. Low coverage is legitimate in real curated
    # records -- Aspergillus nidulans 162425_fgsc-a4_MAT_combined has genes at
    # 3.0%, 18.9% and 30.5% with 100% identity, where the curated CDS span
    # deliberately covers part of a larger deposited protein. Failing those
    # would invalidate correct curation; warning asks the question without
    # overruling a curator.
    if status == "pass" and coverage < WARN_COVERAGE_THRESHOLD:
        status = "warn"
    return MatchScore(percent_identity=round(percent_identity, 2), coverage=round(coverage, 2), status=status)
