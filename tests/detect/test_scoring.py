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
    # FAM_A has 1/2 of its own genes found (g1, attributed to FAM_A).
    # FAM_B has 2/4 of its own genes found (g3, g4, attributed to FAM_B).
    # FAM_B has more genes found in absolute count (2 vs 1) but a lower
    # fraction (0.5 vs 0.5 -- tie) is not useful here, so instead give FAM_B
    # only 1 of its 4 genes found: fraction 0.25 < FAM_A's 0.5. This proves
    # the larger family (more expected genes) does not win on raw hit count;
    # each family is scored strictly against hits attributed to ITS OWN
    # family_key, never against another family's hits.
    cluster = GeneCluster(
        "c1", 1, 100,
        [_hit("g1", FAM_A.key), _hit("g3", FAM_B.key)],
    )
    scores = score_cluster(cluster, [FAM_A, FAM_B])
    by_key = {s.family_key: s for s in scores}
    assert by_key[FAM_A.key].fraction_found == 0.5
    assert by_key[FAM_B.key].fraction_found == 0.25
    assert scores[0].family_key == FAM_A.key  # sorted highest fraction first


def test_score_cluster_does_not_credit_other_family_for_shared_gene_name():
    # FAM_A and FAM_B both declare a gene named "g1". Only FAM_A actually has
    # a hit for it (correctly attributed via family_key=FAM_A.key). FAM_B
    # must NOT be credited for "g1" just because the name string matches --
    # this is the cross-family contamination bug the scoring module exists
    # to prevent.
    cluster = GeneCluster("c1", 1, 100, [_hit("g1", FAM_A.key)])
    scores = score_cluster(cluster, [FAM_A, FAM_B])
    by_key = {s.family_key: s for s in scores}
    assert by_key[FAM_A.key].genes_found == ["g1"]
    # FAM_B has no hits attributed to it at all, so it should not appear
    # in the scores (score_cluster only returns families with >=1 hit).
    assert FAM_B.key not in by_key


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


# Family split across 2 idiomorphs: 2 genes each for MAT1-1/MAT1-2, plus 1
# idiomorph-agnostic flanking gene. Mirrors this session's real-world
# regression (a real single-idiomorph genome scored against the full,
# multi-idiomorph roster was structurally capped below the ambiguity floor).
FAM_IDIOMORPHIC = Family(
    FamilyKey("P", "MAT"), "enum", ["MAT1-1", "MAT1-2"], None,
    [
        {"name": "a1", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-1"]},
        {"name": "a2", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-1"]},
        {"name": "b1", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-2"]},
        {"name": "b2", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-2"]},
        {"name": "flank1", "role": "flanking_conserved"},
    ],
    [1],
)


def test_score_cluster_narrows_to_single_found_idiomorph():
    # Only MAT1-2's genes (b1, b2) plus the idiomorph-agnostic flank1 are
    # found. Against the full 5-gene roster this would be 3/5 = 0.6, but
    # narrowed to MAT1-2 + agnostic genes (b1, b2, flank1) it is 3/3 = 1.0.
    cluster = GeneCluster(
        "c1", 1, 100,
        [_hit("b1", FAM_IDIOMORPHIC.key), _hit("b2", FAM_IDIOMORPHIC.key), _hit("flank1", FAM_IDIOMORPHIC.key)],
    )
    scores = score_cluster(cluster, [FAM_IDIOMORPHIC])
    assert scores[0].fraction_found == 1.0
    assert sorted(scores[0].genes_found) == ["b1", "b2", "flank1"]
    assert scores[0].genes_missing == []


def test_score_cluster_keeps_full_roster_when_both_idiomorphs_found():
    # Genes from BOTH MAT1-1 and MAT1-2 found together -- a real, legitimate
    # homothallic both-idiomorphs-present locus (e.g. curated record
    # db/Ascomycota/Eurotiales/162425_fgsc-a4_MAT_combined, A. nidulans).
    # Must NOT narrow: full 5-gene roster stays in effect.
    cluster = GeneCluster(
        "c1", 1, 100,
        [_hit("a1", FAM_IDIOMORPHIC.key), _hit("b1", FAM_IDIOMORPHIC.key)],
    )
    scores = score_cluster(cluster, [FAM_IDIOMORPHIC])
    assert scores[0].fraction_found == 2 / 5
    assert sorted(scores[0].genes_missing) == ["a2", "b2", "flank1"]


def test_score_cluster_keeps_full_roster_when_no_idiomorph_info_found():
    # Only the idiomorph-agnostic flanking gene is found -- nothing to
    # narrow by, so the full roster is used unchanged.
    cluster = GeneCluster("c1", 1, 100, [_hit("flank1", FAM_IDIOMORPHIC.key)])
    scores = score_cluster(cluster, [FAM_IDIOMORPHIC])
    assert scores[0].fraction_found == 1 / 5
    assert sorted(scores[0].genes_missing) == ["a1", "a2", "b1", "b2"]
