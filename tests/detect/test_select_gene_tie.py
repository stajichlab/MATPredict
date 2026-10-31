"""Two equal alignments of one reference record in a window: the pick does not depend on the tool's output order.

A gene that flanks a MAT locus on both sides gives exonerate two alignments of the same reference protein with the same
score and identity. The tool's order among equal scores is not guaranteed, and the run-to-run result changed with it
(GCA_030583325.1: the sla2 model landed on the left copy in 2 of 10 runs and on the right copy in 8)."""
from __future__ import annotations

import itertools

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.search import _AlignmentRecord, _select_requested_gene

FAMILY = Family(FamilyKey("Ascomycota", "MATsc"), "pattern", None, r"^MAT",
                [{"name": "sla2", "role": "flanking_conserved"}], [4893])
RECORDS = {"recA": FAMILY.key}
ROLES = {FAMILY.key: {"sla2": "flanking_conserved"}}


def _rec(start, strand, query="recA|gene0|sla2", score=500.0):
    feature = ["c1", "exonerate", "gene", str(start), str(start + 266), ".", strand, ".",
               f"gene_id 1 ; sequence {query} ; identity 72.00 ; similarity 80.00"]
    return _AlignmentRecord(query, score, feature, [], ["c1", "exonerate", "similarity", "1", "2", ".", strand, ".",
                                                         "identity 72.00"])


def _chosen(records):
    return _select_requested_gene(records, "sla2", RECORDS, ROLES)[0].feature[3]


def test_equal_copies_resolve_to_the_lower_start_in_either_order():
    left, right = _rec(40124, "+"), _rec(56610, "-")
    assert _chosen([left, right]) == "40124"
    assert _chosen([right, left]) == "40124"


def test_equal_copies_at_one_start_resolve_by_strand_in_any_order():
    a, b = _rec(100, "+"), _rec(100, "-")
    picks = {_chosen(list(p)) + _select_requested_gene(list(p), "sla2", RECORDS, ROLES)[0].feature[6]
             for p in itertools.permutations([a, b])}
    assert len(picks) == 1


def test_a_higher_score_still_wins_over_position():
    assert _chosen([_rec(40124, "+", score=400.0), _rec(56610, "-", score=500.0)]) == "56610"
