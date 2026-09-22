"""`idiomorph_margin` must describe the margin of the call that was actually made.

Curator's ruling, 2026-09-21: when a `partial_locus` resolves to an idiomorph,
"report the margin".

It was already reported -- but it described the WRONG mechanism. Measured on
Actinomucor sp. NRRL A-23671 from the ground-truth set:

    idiomorph            Plus
    idiomorph_margin     1.077                       <- identity, from the
                                                        sexM/sexP OVERLAP
                                                        resolution
    idiomorph_candidates Plus 35.8, Minus 31.2       <- bitscore, from the VOTE
                                                        that actually decided

Since the 2026-09-21 change the idiomorph is decided by the vote, not by the
overlap resolution, so a field named `idiomorph_margin` reporting 1.077 while
the decision turned on a gap of 4.6 is actively misleading -- a reader
filtering on a thin margin would filter the wrong thing.

The margin now comes from the same ranking that made the call. The overlap
resolution's own margin is unchanged and still visible in
`idiomorph_resolutions`, where it belongs.
"""
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.idiomorph import idiomorph_candidates, idiomorph_margin_from_vote
from MATPredict.detect.search import SearchHit

KEY = FamilyKey("Mucoromycota", "MAT")
FAM = Family(
    key=KEY, vocabulary_type="enum", idiomorph_values=["Plus", "Minus"],
    idiomorph_pattern=None,
    genes=[
        {"name": "sexP", "role": "core_MAT", "present_in_idiomorphs": ["Plus"]},
        {"name": "sexM", "role": "core_MAT", "present_in_idiomorphs": ["Minus"]},
    ],
    taxonomic_scope=[4827],
)


def _h(gene, bits):
    return SearchHit(KEY, gene, "core_MAT", "c1", 100, 400, "+", 50.0, "r",
                     "tblastn_genome", bitscore=bits)


def test_the_margin_is_the_gap_between_the_top_two():
    hits = [_h("sexP", 35.8), _h("sexM", 31.2)]
    ranked = idiomorph_candidates(FAM, ["sexP", "sexM"], hits)
    assert idiomorph_margin_from_vote(ranked) == 4.6


def test_a_single_candidate_has_no_margin():
    """Nothing to be narrow against. None, not zero -- zero would read as a tie."""
    ranked = idiomorph_candidates(FAM, ["sexP"], [_h("sexP", 35.8)])
    assert idiomorph_margin_from_vote(ranked) is None


def test_no_candidates_has_no_margin():
    assert idiomorph_margin_from_vote([]) is None


def test_an_exact_tie_is_a_zero_margin():
    hits = [_h("sexP", 40.0), _h("sexM", 40.0)]
    ranked = idiomorph_candidates(FAM, ["sexP", "sexM"], hits)
    assert idiomorph_margin_from_vote(ranked) == 0.0
