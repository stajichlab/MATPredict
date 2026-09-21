# tests/detect/test_optional_genes.py
"""An optional gene is reported when found and never counted against a locus.

Curator ruling, J. Stajich, 2026-09-20: "btbA shouldn't count in the
denominator but we should know when it is present."

The measurement behind it, from the 283-genome BFD Mucoromycota sweep:
`btbA` is present in 0 of 207 Minus loci and only 57 of 220 Plus loci (26%).
It is perfectly specific to Plus and weakly sensitive -- presence argues
Plus, absence argues nothing.

Because `btbA` carries `present_in_idiomorphs: ["Plus"]`,
`expected_genes_for_idiomorph` put it in the expected roster for every Plus
locus. A Plus locus was therefore scored against 6 expected genes and a Minus
locus against 5, while 74% of Plus loci lack `btbA` -- an asymmetric
`fraction_found` penalty that no Minus locus could incur.

This is the same idea as `genes_not_searchable`, which already exists: a gene
that cannot count against the locus is dropped from the denominator rather
than reported as missing. The difference is the reason -- unsearchable means
"the run could not have found it", optional means "the biology does not
require it".
"""
from __future__ import annotations

from pathlib import Path

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import (
    Family,
    FamilyKey,
    load_all_families,
)
from MATPredict.detect.scoring import score_cluster
from MATPredict.detect.search import SearchHit

KEY = FamilyKey("Mucoromycota", "MAT")
FAM = Family(
    KEY, "enum", ["Plus", "Minus"], None,
    [
        {"name": "tptA", "role": "flanking_conserved"},
        {"name": "sexP", "role": "core_MAT", "present_in_idiomorphs": ["Plus"]},
        {"name": "sexM", "role": "core_MAT", "present_in_idiomorphs": ["Minus"]},
        {"name": "rnhA", "role": "flanking_conserved"},
        {"name": "algA", "role": "flanking_variable"},
        {"name": "glrA", "role": "flanking_variable"},
        {"name": "btbA", "role": "flanking_variable",
         "present_in_idiomorphs": ["Plus"], "optional": True},
    ],
    [4827],
)
_ROLES = {g["name"]: g["role"] for g in FAM.genes}
ALL = {g["name"] for g in FAM.genes}


def _cluster(*genes):
    hits = [
        SearchHit(KEY, g, _ROLES[g], "c1", 1000 + 5000 * i, 3000 + 5000 * i,
                  "+", 60.0, "rec1", "diamond_proteome")
        for i, g in enumerate(genes)
    ]
    return GeneCluster("c1", hits[0].start, hits[-1].end, hits)


def _score(*genes):
    return score_cluster(_cluster(*genes), [FAM], searchable_genes={KEY: ALL})[0]


def test_an_absent_optional_gene_does_not_reduce_the_fraction():
    """The whole point. A complete Plus locus with no btbA must score 1.0,
    not 5/6 = 0.833."""
    s = _score("tptA", "sexP", "rnhA", "algA", "glrA")
    assert s.fraction_found == 1.0


def test_an_absent_optional_gene_is_not_reported_missing():
    s = _score("tptA", "sexP", "rnhA", "algA", "glrA")
    assert "btbA" not in s.genes_missing


def test_a_present_optional_gene_is_reported():
    """"We should know when it is present" -- it stays in genes_found and is
    named separately so the arithmetic stays auditable."""
    s = _score("tptA", "sexP", "rnhA", "algA", "glrA", "btbA")
    assert "btbA" in s.genes_found
    assert s.genes_optional_found == ["btbA"]


def test_a_present_optional_gene_cannot_push_the_fraction_over_one():
    s = _score("tptA", "sexP", "rnhA", "algA", "glrA", "btbA")
    assert s.fraction_found == 1.0


def test_an_optional_gene_does_not_rescue_an_otherwise_poor_locus():
    """btbA must not substitute for a real gene. sexP + btbA is one required
    gene of five, whatever btbA adds."""
    s = _score("sexP", "btbA")
    assert s.fraction_found == 1 / 5
    assert s.genes_optional_found == ["btbA"]


def test_plus_and_minus_now_share_a_denominator_size():
    """The asymmetry this ruling removes: before, Plus was scored against 6
    expected genes and Minus against 5."""
    plus = _score("sexP")
    minus = _score("sexM")
    assert len(plus.genes_found) + len(plus.genes_missing) == 5
    assert len(minus.genes_found) + len(minus.genes_missing) == 5


def test_the_curated_mucoromycota_btba_is_marked_optional():
    """The ruling has to reach the real curation data, not just the schema."""
    db = Path(__file__).resolve().parents[2] / "db"
    fam = [f for f in load_all_families(db) if f.key == KEY][0]
    btba = [g for g in fam.genes if g["name"] == "btbA"][0]
    assert btba.get("optional") is True
    # And nothing else is optional without a ruling behind it.
    assert [g["name"] for g in fam.genes if g.get("optional")] == ["btbA"]
