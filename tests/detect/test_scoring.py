# tests/detect/test_scoring.py
from __future__ import annotations

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.search import SearchHit
from MATPredict.detect.scoring import score_cluster, is_ambiguous

FAM_A = Family(FamilyKey("P", "A"), "enum", ["a", "alpha"], None,
               [{"name": "g1", "role": "core_MAT"}, {"name": "g2", "role": "core_MAT"}], [1])
FAM_B = Family(FamilyKey("P", "B"), "enum", ["a", "alpha"], None,
               [{"name": "g1", "role": "core_MAT"}, {"name": "g3", "role": "core_MAT"},
                {"name": "g4", "role": "core_MAT"}, {"name": "g5", "role": "core_MAT"}], [1])


def _hit(gene_name, family_key):
    return SearchHit(family_key, gene_name, "core_MAT", "c1", 1, 100, "+", 90.0, "rec1", "diamond_proteome")


def test_score_cluster_is_fractional_not_summed():
    # FAM_A: 1/2 genes found. FAM_B: 1/4 genes found. Gene count alone must not win for FAM_B.
    cluster = GeneCluster("c1", 1, 100, [_hit("g1", FAM_A.key)])
    scores = score_cluster(cluster, [FAM_A, FAM_B])
    by_key = {s.family_key: s for s in scores}
    assert by_key[FAM_A.key].fraction_found == 0.5
    assert by_key[FAM_B.key].fraction_found == 0.25
    assert scores[0].family_key == FAM_A.key  # sorted highest fraction first


def test_score_cluster_reports_missing_genes():
    cluster = GeneCluster("c1", 1, 100, [_hit("g1", FAM_A.key)])
    scores = score_cluster(cluster, [FAM_A])
    assert scores[0].genes_found == ["g1"]
    assert scores[0].genes_missing == ["g2"]


def test_is_ambiguous_detects_multiple_high_scorers():
    from MATPredict.detect.scoring import FamilyScore
    scores = [
        FamilyScore(FAM_A.key, 0.9, ["g1", "g2"], []),
        FamilyScore(FAM_B.key, 0.75, ["g1", "g3", "g4"], ["g5"]),
    ]
    assert is_ambiguous(scores, floor=0.5) is True


def test_is_ambiguous_false_with_one_clear_winner():
    from MATPredict.detect.scoring import FamilyScore
    scores = [FamilyScore(FAM_A.key, 0.9, ["g1", "g2"], []), FamilyScore(FAM_B.key, 0.1, ["g1"], ["g3", "g4", "g5"])]
    assert is_ambiguous(scores, floor=0.5) is False
