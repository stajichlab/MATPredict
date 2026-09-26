# tests/detect/test_scoring_slot.py
"""A `slot` groups idiomorph-alternative optional genes (SXI1/SXI2) into ONE
bonus scoring unit: it enters numerator and denominator only when a member is
found, so it can raise a score but never lower one.

Measured case (results/2026-09-26_cryptococcus_zero/NOTE.md): four MATalpha
C. neoformans genomes have the MAT locus split over 5-6 contigs. The true piece
holds SXI1 + FAO1. With SXI1 merely optional it scored 1/3 (FAO1 of FAO1, PAN6,
STE3) and fell below the 0.5 floor.
"""
from __future__ import annotations

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.scoring import score_cluster
from MATPredict.detect.search import SearchHit

KEY = FamilyKey("Basidiomycota", "MAT")
GENES = [
    {"name": "SXI1", "role": "core_MAT", "present_in_idiomorphs": ["alpha"],
     "optional": True, "slot": "SXI"},
    {"name": "SXI2", "role": "core_MAT", "present_in_idiomorphs": ["a"],
     "optional": True, "slot": "SXI"},
    {"name": "RPL39", "role": "core_MAT", "optional": True},
    {"name": "FAO1", "role": "flanking_conserved"},
    {"name": "PAN6", "role": "flanking_conserved"},
    {"name": "STE3", "role": "core_MAT"},
]
FAM = Family(KEY, "enum", ["alpha", "a"], None, GENES, [5234])
NO_SLOT = Family(KEY, "enum", ["alpha", "a"], None,
                 [{k: v for k, v in g.items() if k != "slot"} for g in GENES], [5234])


def _hit(gene):
    return SearchHit(KEY, gene, "core_MAT", "c1", 1, 100, "+", 90.0, "rec1", "tblastn")


def _score(fam, genes):
    cluster = GeneCluster("c1", 1, 100, [_hit(g) for g in genes])
    (score,) = score_cluster(cluster, [fam])
    return score


def test_found_slot_member_counts_once_in_numerator_and_denominator():
    s = _score(FAM, ["SXI1", "FAO1"])
    assert s.fraction_found == 2 / 4
    assert "SXI1" in s.genes_found
    assert "SXI1" not in s.genes_optional_found


def test_without_slot_the_same_cluster_scores_one_third():
    assert _score(NO_SLOT, ["SXI1", "FAO1"]).fraction_found == 1 / 3


def test_absent_slot_never_lowers_a_score():
    for genes in (["FAO1"], ["FAO1", "STE3"], ["FAO1", "PAN6", "STE3"]):
        assert _score(FAM, genes).fraction_found == _score(NO_SLOT, genes).fraction_found


def test_two_slot_members_in_one_cluster_add_nothing():
    # SXI1 and SXI2 are alternatives: both hitting ONE cluster is what a
    # homeodomain paralog does, not evidence of either idiomorph. Measured on
    # 243 Cryptococcus genomes (results/2026-09-26_sxi_slot_gateA/): 12 new
    # calls held both, 3 of them turning a single-idiomorph genome into a+alpha
    # on another chromosome. So the slot scores nothing then -- the score is
    # exactly what it was with plain optional genes.
    s = _score(FAM, ["SXI1", "SXI2", "FAO1"])
    assert s.fraction_found == _score(NO_SLOT, ["SXI1", "SXI2", "FAO1"]).fraction_found
    assert "SXI1" in s.genes_found and "SXI2" in s.genes_found


def test_slot_alone_cannot_manufacture_a_score():
    only_optional = Family(KEY, "enum", ["alpha", "a"], None,
                           [g for g in GENES if g.get("optional")], [5234])
    assert _score(only_optional, ["SXI1"]).fraction_found == 0.0


def test_plain_optional_gene_is_still_unscored():
    s = _score(FAM, ["RPL39", "FAO1"])
    assert s.fraction_found == 1 / 3
    assert s.genes_optional_found == ["RPL39"]
