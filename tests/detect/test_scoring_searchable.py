# tests/detect/test_scoring_searchable.py
"""The `fraction_found` denominator counts only genes this run could have found.

A gene named in a locus's `order.yml` roster but absent from the reference set
the run actually searched with can never be found, by any genome. Counting it
in the denominator does not make the score conservative, it makes it
uninterpretable: the maximum achievable value drifts with gaps in the curated
database rather than with the biology being measured.

This is not hypothetical. Measured on the live database, Mucoromycota's MAT
roster names 7 genes, of which `algA` and `glrA` have zero reference proteins
anywhere in `db/`. Every Plus genome in the 23-genome ground-truth set can
therefore reach at most 4/6 once its idiomorph roster is narrowed -- and after
overlapping sexM/sexP hits are resolved to one gene it lands on exactly 0.500,
the rejection boundary, surviving only because the floor test is `<` and not
`<=`.
"""
from __future__ import annotations

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.scoring import score_cluster
from MATPredict.detect.search import SearchHit

FAM = Family(
    FamilyKey("P", "A"), "enum", ["Plus", "Minus"], None,
    [{"name": "g1", "role": "core_MAT"}, {"name": "g2", "role": "flanking_variable"}],
    [1],
)


def _hit(gene_name, family_key, role="core_MAT"):
    return SearchHit(
        family_key, gene_name, role, "c1", 1, 100, "+", 90.0, "rec1", "diamond_proteome"
    )


def test_denominator_excludes_genes_with_no_reference_protein():
    # g2 is in the roster but has no reference protein in the searched set, so
    # a cluster holding every findable gene scores 1.0, not 0.5.
    cluster = GeneCluster("c1", 1, 100, [_hit("g1", FAM.key)])
    scores = score_cluster(cluster, [FAM], searchable_genes={FAM.key: {"g1"}})
    assert scores[0].fraction_found == 1.0
    assert scores[0].genes_found == ["g1"]
    assert scores[0].genes_missing == []
    assert scores[0].genes_not_searchable == ["g2"]


def test_a_searchable_gene_that_was_not_found_still_counts_as_missing():
    # The exclusion must apply ONLY to genes that could not be searched for.
    # A gene with a reference protein that simply was not found is a real
    # absence and must keep lowering the score, or the denominator becomes a
    # way to make any locus look complete.
    cluster = GeneCluster("c1", 1, 100, [_hit("g1", FAM.key)])
    scores = score_cluster(cluster, [FAM], searchable_genes={FAM.key: {"g1", "g2"}})
    assert scores[0].fraction_found == 0.5
    assert scores[0].genes_missing == ["g2"]
    assert scores[0].genes_not_searchable == []


def test_omitting_searchable_genes_keeps_the_whole_roster_searchable():
    # Every existing caller passes nothing, and must keep its exact behaviour:
    # absence of the argument means "no information about searchability", not
    # "nothing is searchable".
    cluster = GeneCluster("c1", 1, 100, [_hit("g1", FAM.key)])
    scores = score_cluster(cluster, [FAM])
    assert scores[0].fraction_found == 0.5
    assert scores[0].genes_missing == ["g2"]
    assert scores[0].genes_not_searchable == []


def test_a_family_absent_from_the_map_keeps_its_whole_roster():
    # A map that simply has no entry for this family carries no information
    # about it. Treating a missing key as "nothing searchable" would divide by
    # zero; treating it as "all searchable" preserves the prior behaviour.
    cluster = GeneCluster("c1", 1, 100, [_hit("g1", FAM.key)])
    scores = score_cluster(cluster, [FAM], searchable_genes={FamilyKey("P", "other"): {"x"}})
    assert scores[0].fraction_found == 0.5
    assert scores[0].genes_not_searchable == []


def test_searchability_is_applied_after_the_idiomorph_roster_is_narrowed():
    # Narrowing to one idiomorph and dropping unsearchable genes compose: the
    # denominator is the intersection of "this idiomorph expects it" and "we
    # could have found it". Getting the order wrong would report the other
    # idiomorph's genes as not searchable.
    family = Family(
        FamilyKey("P", "B"), "enum", ["Plus", "Minus"], None,
        [
            {"name": "flank", "role": "flanking_conserved"},
            {"name": "p_only", "role": "core_MAT", "present_in_idiomorphs": ["Plus"]},
            {"name": "m_only", "role": "core_MAT", "present_in_idiomorphs": ["Minus"]},
            {"name": "p_extra", "role": "flanking_variable", "present_in_idiomorphs": ["Plus"]},
        ],
        [1],
    )
    cluster = GeneCluster(
        "c1", 1, 100,
        [_hit("flank", family.key, "flanking_conserved"), _hit("p_only", family.key)],
    )
    # p_extra is Plus-expected but has no reference protein.
    scores = score_cluster(
        cluster, [family], searchable_genes={family.key: {"flank", "p_only", "m_only"}}
    )
    assert scores[0].fraction_found == 1.0
    assert scores[0].genes_found == ["flank", "p_only"]
    assert scores[0].genes_not_searchable == ["p_extra"]
    # m_only belongs to the other idiomorph; it is neither missing nor
    # unsearchable here, it is simply not expected.
    assert "m_only" not in scores[0].genes_missing
    assert "m_only" not in scores[0].genes_not_searchable
