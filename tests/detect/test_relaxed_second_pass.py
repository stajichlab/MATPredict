# tests/detect/test_relaxed_second_pass.py
"""A second, lower-bar pass, run only when the strict pass finds nothing.

Curator's ruling, 2026-09-19/20: a genuinely fragmented locus -- a real MAT
locus split by a contig break, so its flanking genes sit on other contigs --
can fall below the strict fraction floor through no fault of its own. The
relaxed pass admits it on "a core gene plus another gene, not necessarily a
flank", which is character-for-character the existing `EvidenceFloor`
(min_hits=2, require_core_role=True) already used to admit clusters to
polishing. One curator-ruled bar serving two purposes, so there is only one
number to defend.

It runs ONLY when the strict pass found nothing genome-wide, so a genome with
any confident call is completely unaffected.

Scope note, measured: for Mucoromycota this admits nothing new. With the
searchable-only denominators (Plus 4, Minus 3), "core + one other gene"
already clears the strict 0.50 floor -- 2/4 = 0.50 and 2/3 = 0.667. The pass
earns its place in phyla with richer rosters, where a real core+flank pair can
score 2/8 = 0.25 and be rejected. That is why these tests use an 8-gene family.
"""
from __future__ import annotations

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.pipeline import EvidenceFloor, _relaxed_results
from MATPredict.detect.scoring import score_cluster
from MATPredict.detect.search import SearchHit

RICH = Family(
    FamilyKey("Ascomycota", "MAT"), "enum", ["MAT1-1", "MAT1-2"], None,
    [{"name": f"g{i}", "role": "core_MAT" if i == 1 else "flanking_variable"}
     for i in range(1, 9)],
    [4890],
)
SEARCHABLE = {RICH.key: {f"g{i}" for i in range(1, 9)}}


def _hit(gene, role, identity=60.0, start=1, end=100):
    return SearchHit(
        RICH.key, gene, role, "c1", start, end, "+", identity, "rec1",
        "diamond_proteome",
    )


def _cluster(*hits):
    return GeneCluster("c1", 1, 1000, list(hits))


def _core_plus_flank():
    return _cluster(
        _hit("g1", "core_MAT", 92.6, 1, 100),
        _hit("g2", "flanking_variable", 88.0, 400, 500),
    )


def test_a_core_plus_one_other_gene_is_admitted_when_nothing_else_was_found():
    # 2 of 8 expected = 0.25, well under the strict 0.50 floor, yet it is a
    # real core gene and a real flank at high identity.
    cluster = _core_plus_flank()
    results = _relaxed_results(
        [cluster], [RICH], searchable_genes=SEARCHABLE,
        evidence_floor=EvidenceFloor(),
    )
    assert [r.family_key for r in results] == [RICH.key]


def test_the_relaxed_call_is_labelled_and_capped():
    results = _relaxed_results(
        [_core_plus_flank()], [RICH], searchable_genes=SEARCHABLE,
        evidence_floor=EvidenceFloor(),
    )
    assert results[0].detection_pass == "relaxed"
    assert results[0].confidence in ("medium", "low")


def test_a_lone_core_gene_is_still_rejected():
    # The bar is a COUNT of distinct genes, so one gene never passes however
    # good it looks. This is what keeps the 22 lone-HMG-gene genera out.
    cluster = _cluster(_hit("g1", "core_MAT", 99.0))
    assert _relaxed_results(
        [cluster], [RICH], searchable_genes=SEARCHABLE,
        evidence_floor=EvidenceFloor(),
    ) == []


def test_two_genes_with_no_core_are_still_rejected():
    cluster = _cluster(
        _hit("g2", "flanking_variable", 90.0, 1, 100),
        _hit("g3", "flanking_variable", 88.0, 400, 500),
    )
    assert _relaxed_results(
        [cluster], [RICH], searchable_genes=SEARCHABLE,
        evidence_floor=EvidenceFloor(),
    ) == []


def test_a_cluster_that_already_clears_the_strict_floor_is_not_re_reported():
    # The relaxed pass only ever runs when the strict pass produced nothing,
    # so anything it reports is by construction sub-floor. A cluster scoring
    # at or above the floor would be a duplicate of a strict result.
    cluster = _cluster(*[
        _hit(f"g{i}", "core_MAT" if i == 1 else "flanking_variable", 80.0, i * 100, i * 100 + 50)
        for i in range(1, 6)
    ])
    assert score_cluster(cluster, [RICH], searchable_genes=SEARCHABLE)[0].fraction_found >= 0.5
    assert _relaxed_results(
        [cluster], [RICH], searchable_genes=SEARCHABLE,
        evidence_floor=EvidenceFloor(),
    ) == []


def test_a_strict_result_reports_detection_pass_strict():
    # The field must be present on every result, not only relaxed ones, or a
    # reader cannot tell "strict" from "this build predates the field".
    from MATPredict.detect.pipeline import DetectionResult

    r = DetectionResult(
        family_key=RICH.key, contig="c1", start=1, end=100, confidence="high",
        idiomorph="MAT1-1", ambiguous_with=[], genes_found=["g1"],
        genes_missing=[], fragmented=False,
    )
    assert r.detection_pass == "strict"
