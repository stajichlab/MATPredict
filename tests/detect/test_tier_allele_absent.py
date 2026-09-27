# tests/detect/test_tier_allele_absent.py
"""Allele-absent genes and the confidence tier (curator's ruling 2026-09-27,
variant B' with the closeness guard; replay in
results/2026-09-27_tier_rule_replay).

A gene whose `present_in_idiomorphs` excludes the called idiomorph is ignored
by the tier -- in the unpolished check and in the found set that decides the
expected core -- only when it is weak (< 50% identity) AND at least 10 points
below the called allele's best modelled core gene."""
from __future__ import annotations

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.scoring import FamilyScore
from MATPredict.detect.tiering import (
    ALLELE_ABSENT_MAX_IDENTITY,
    ALLELE_ABSENT_MIN_GAP,
    allele_absent_genes_to_ignore,
    assign_tier,
)

# Shaped on the Wallemia wallMAT roster: v1 carries SXI1/HMG/STE3, v2 carries
# only STE3v2; BAP31/CAF1 flank both and are optional.
WALL = Family(
    FamilyKey("Basidiomycota", "wallMAT"), "enum", ["v1", "v2"], None,
    [
        {"name": "SXI1", "role": "core_MAT", "present_in_idiomorphs": ["v1"]},
        {"name": "HMG", "role": "core_MAT", "present_in_idiomorphs": ["v1"]},
        {"name": "STE3", "role": "core_MAT", "present_in_idiomorphs": ["v1"]},
        {"name": "STE3v2", "role": "core_MAT", "present_in_idiomorphs": ["v2"]},
        {"name": "BAP31", "role": "flanking_variable", "optional": True},
        {"name": "CAF1", "role": "flanking_variable", "optional": True},
    ],
    [431958],
)


def test_thresholds_are_the_ruled_values():
    assert ALLELE_ABSENT_MAX_IDENTITY == 50.0
    assert ALLELE_ABSENT_MIN_GAP == 10.0


def test_weak_distant_cross_hit_is_ignored():
    ignore = allele_absent_genes_to_ignore(
        WALL, "v2", {"STE3v2": 88.0, "HMG": 36.0, "BAP31": 90.0}, modelled={"STE3v2", "BAP31"},
    )
    assert ignore == frozenset({"HMG"})


def test_cross_hit_at_or_above_50_is_kept():
    ignore = allele_absent_genes_to_ignore(
        WALL, "v2", {"STE3v2": 95.0, "HMG": 60.5}, modelled={"STE3v2"},
    )
    assert ignore == frozenset()


def test_cross_hit_close_to_the_called_gene_is_kept():
    # 44 vs 46: the "absent" allele is nearly as strong as the called one.
    ignore = allele_absent_genes_to_ignore(
        WALL, "v2", {"STE3v2": 46.0, "HMG": 44.0}, modelled={"STE3v2"},
    )
    assert ignore == frozenset()


def test_gap_of_exactly_ten_points_is_ignored():
    ignore = allele_absent_genes_to_ignore(
        WALL, "v2", {"STE3v2": 45.0, "HMG": 35.0}, modelled={"STE3v2"},
    )
    assert ignore == frozenset({"HMG"})


def test_undetermined_call_ignores_nothing():
    assert allele_absent_genes_to_ignore(
        WALL, "undetermined", {"STE3v2": 88.0, "HMG": 36.0}, modelled={"STE3v2"},
    ) == frozenset()
    assert allele_absent_genes_to_ignore(
        WALL, None, {"STE3v2": 88.0, "HMG": 36.0}, modelled={"STE3v2"},
    ) == frozenset()


def test_no_modelled_core_of_the_called_allele_ignores_nothing():
    # The guard compares against a MODELLED gene of the called allele; with
    # none there is nothing to measure the gap from.
    assert allele_absent_genes_to_ignore(
        WALL, "v2", {"STE3v2": 88.0, "HMG": 36.0}, modelled={"BAP31"},
    ) == frozenset()


def test_gene_without_identity_is_kept():
    assert allele_absent_genes_to_ignore(
        WALL, "v2", {"STE3v2": 88.0, "HMG": None}, modelled={"STE3v2"},
    ) == frozenset()


def test_one_strong_other_allele_gene_blocks_ignoring_any():
    # Guard applied per CALL: a weak HMG cross-hit is not ignored while the
    # call also carries a strong SXI1 (the collapsed-locus shape).
    assert allele_absent_genes_to_ignore(
        WALL, "v2", {"STE3v2": 88.0, "HMG": 30.0, "SXI1": 55.0}, modelled={"STE3v2"},
    ) == frozenset()


def test_one_close_other_allele_gene_blocks_ignoring_any():
    # Serinales shape: called best 41.2, weak MTLA2-like 28.8 but a kept
    # 45.3 other-allele gene -- nothing is ignored.
    assert allele_absent_genes_to_ignore(
        WALL, "v2", {"STE3v2": 41.2, "HMG": 28.8, "SXI1": 45.3}, modelled={"STE3v2"},
    ) == frozenset()


def test_all_weak_distant_other_allele_genes_are_ignored_together():
    assert allele_absent_genes_to_ignore(
        WALL, "v2", {"STE3v2": 88.0, "HMG": 36.0, "STE3": 40.0}, modelled={"STE3v2"},
    ) == frozenset({"HMG", "STE3"})


def test_agnostic_genes_are_never_ignored():
    # BAP31 has no present_in_idiomorphs: it belongs to every allele.
    ignore = allele_absent_genes_to_ignore(
        WALL, "v2", {"STE3v2": 88.0, "BAP31": 30.0}, modelled={"STE3v2"},
    )
    assert ignore == frozenset()


def test_ignored_cross_hit_no_longer_widens_the_expected_core():
    # The Wallemia v2 case: a weak HMG cross-hit makes the found genes name
    # BOTH alleles, the expected core becomes the full roster, SXI1 is
    # "missing" and the call is capped at medium.
    score = FamilyScore(WALL.key, 1.0, ["STE3v2", "HMG", "BAP31", "CAF1"], [])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, WALL, cluster, any_gene_unpolished=False, fragmented=False) == "medium"
    assert assign_tier(
        score, WALL, cluster, any_gene_unpolished=False, fragmented=False,
        ignore_genes=frozenset({"HMG"}),
    ) == "high"


def test_ignore_genes_defaults_to_current_behaviour():
    score = FamilyScore(WALL.key, 1.0, ["STE3v2", "BAP31"], [])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, WALL, cluster, any_gene_unpolished=False, fragmented=False) == "high"


def test_homothallic_candidate_ignores_nothing():
    # Curator's ruling 2026-09-27: a homothallic call rests on BOTH alleles by
    # definition, so neither may be dropped from the tier. Case that forced it:
    # Serinales GCA_030462985.1, homothallic_candidate, raised to high by
    # ignoring its MTLA2 (38.75%).
    args = (WALL, "v2", {"STE3v2": 88.0, "HMG": 36.0}, {"STE3v2"})
    assert allele_absent_genes_to_ignore(*args) == frozenset({"HMG"})
    assert allele_absent_genes_to_ignore(
        *args, locus_class="homothallic_candidate",
    ) == frozenset()
    assert allele_absent_genes_to_ignore(
        *args, locus_class="mat_locus",
    ) == frozenset({"HMG"})
