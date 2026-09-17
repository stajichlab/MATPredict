# tests/detect/test_tiering.py
from __future__ import annotations

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.scoring import FamilyScore
from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.tiering import assign_tier, has_flanking_conserved

FLANKLESS = Family(FamilyKey("P", "aLocus"), "pattern", None, "^a[0-9]+$",
                    [{"name": "mfa1", "role": "core_MAT"}, {"name": "pra1", "role": "core_MAT"}], [1])
FLANKED = Family(FamilyKey("P", "MAT"), "enum", ["a", "alpha"], None,
                  [{"name": "STE3", "role": "core_MAT"}, {"name": "flank1", "role": "flanking_conserved"}], [1])


def test_has_flanking_conserved():
    assert has_flanking_conserved(FLANKED) is True
    assert has_flanking_conserved(FLANKLESS) is False


def test_flankless_family_high_tier_on_all_core_genes_found():
    score = FamilyScore(FLANKLESS.key, 1.0, ["mfa1", "pra1"], [])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, FLANKLESS, cluster, second_pass_used=False, fragmented=False) == "high"


def test_flanked_family_medium_tier_when_core_gene_only_found_via_second_pass():
    score = FamilyScore(FLANKED.key, 1.0, ["STE3", "flank1"], [])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, FLANKED, cluster, second_pass_used=True, fragmented=False) == "medium"


def test_flanked_family_medium_tier_when_no_flanking_gene_found():
    score = FamilyScore(FLANKED.key, 0.5, ["STE3"], ["flank1"])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, FLANKED, cluster, second_pass_used=False, fragmented=False) == "medium"


def test_partial_match_is_medium():
    score = FamilyScore(FLANKLESS.key, 0.5, ["mfa1"], ["pra1"])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, FLANKLESS, cluster, second_pass_used=False, fragmented=False) == "medium"


def test_fragmented_locus_downgraded_one_tier():
    score = FamilyScore(FLANKLESS.key, 1.0, ["mfa1", "pra1"], [])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, FLANKLESS, cluster, second_pass_used=False, fragmented=True) == "medium"
