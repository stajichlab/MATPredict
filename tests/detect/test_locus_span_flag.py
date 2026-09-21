# tests/detect/test_locus_span_flag.py
"""An implausibly wide locus is FLAGGED, never dropped.

Curator ruling, J. Stajich, 2026-09-20, after the 283-genome BFD Mucoromycota
sweep reported spans up to 154,355 bp against ground-truth loci of
6,795-13,089 bp: large loci are real and must still be found -- this is
consistent with unpublished findings by a former graduate student -- so the
bound is a REPORTING flag, not a filter. Set at 200 kb as the widest a
Mucoromycota MAT locus is currently expected to be.

Measured consequence on that sweep: 0 of 687 loci exceed 200 kb (max 154,355),
so this flag fires on nothing already observed. It is a guard-rail against a
runaway cluster, not a filter on present output. Had it been set at the
120 kb first discussed, it would have flagged 3 loci -- including a
Blakeslea trispora call with six genes at 87.2% identity, which is very
likely real and sits only 473 bp over that line.
"""
from __future__ import annotations

from MATPredict.detect.family_registry import (
    DEFAULT_MAX_PLAUSIBLE_LOCUS_SPAN_BP,
    Family,
    FamilyKey,
)

FAM = Family(
    FamilyKey("Mucoromycota", "MAT"), "enum", ["Plus", "Minus"], None,
    [
        {"name": "sexP", "role": "core_MAT", "present_in_idiomorphs": ["Plus"]},
        {"name": "rnhA", "role": "flanking_conserved"},
    ],
    [4827],
)


def test_the_default_bound_is_200kb():
    assert DEFAULT_MAX_PLAUSIBLE_LOCUS_SPAN_BP == 200_000


def test_a_family_carries_the_bound_as_curation_data():
    """On the family, like max_cluster_gap_bp -- how wide a MAT locus gets is
    a property of the clade's architecture, not a global tuning knob."""
    assert FAM.max_plausible_locus_span_bp == 200_000
    wide = Family(
        FAM.key, "enum", ["Plus", "Minus"], None, FAM.genes, [4827],
        max_plausible_locus_span_bp=500_000,
    )
    assert wide.max_plausible_locus_span_bp == 500_000


def test_an_ordinary_locus_is_not_flagged():
    from MATPredict.detect.pipeline import span_exceeds_plausible_bound

    # The widest locus actually seen on the 283-genome sweep.
    assert span_exceeds_plausible_bound(1, 154_355, FAM) is False


def test_a_locus_wider_than_the_bound_is_flagged():
    from MATPredict.detect.pipeline import span_exceeds_plausible_bound

    assert span_exceeds_plausible_bound(1, 200_002, FAM) is True


def test_the_flag_reaches_the_report_and_defaults_to_false():
    """A flag nobody can see is not a flag. It must be on every result, so a
    reader can tell "not flagged" from "this build predates the field"."""
    from MATPredict.detect.pipeline import DetectionResult
    from MATPredict.detect.report import _result_doc

    r = DetectionResult(
        family_key=FAM.key, contig="c1", start=1, end=100, confidence="high",
        idiomorph="Plus", ambiguous_with=[], genes_found=["sexP"],
        genes_missing=[], fragmented=False,
    )
    assert r.span_exceeds_plausible_bound is False
    assert _result_doc(r)["span_exceeds_plausible_bound"] is False

    wide = DetectionResult(
        family_key=FAM.key, contig="c1", start=1, end=300_000,
        confidence="high", idiomorph="Plus", ambiguous_with=[],
        genes_found=["sexP"], genes_missing=[], fragmented=False,
        span_exceeds_plausible_bound=True,
    )
    assert _result_doc(wide)["span_exceeds_plausible_bound"] is True


def test_the_bound_is_inclusive_so_exactly_at_it_is_not_flagged():
    """A locus exactly at the bound is within what the curator called
    acceptable, so `>` and not `>=`. The Blakeslea case is why this edge
    matters: a call 473 bp over a threshold should not be treated
    differently from one 473 bp under it without a reason."""
    from MATPredict.detect.pipeline import span_exceeds_plausible_bound

    assert span_exceeds_plausible_bound(1, 200_000, FAM) is False
