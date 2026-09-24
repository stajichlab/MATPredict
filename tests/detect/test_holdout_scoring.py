"""Scoring one held-out record: found, right idiomorph, withheld, or unfindable.

The first scorer counted ANY overlapping locus as a hit and never looked at the
idiomorph, so S. cerevisiae HMRa scored a hit at record radius while the locus
there was called MATalpha. It also had no way to separate the two kinds of
"no reference" the doc tables relied on; those came from a manual step nobody
could reproduce. And a locus the modelled-gene bar withheld at the right place
looked the same as no call at all.
"""
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.holdout import score_holdout

MAT = FamilyKey("Ascomycota", "MAT")
MATSC = FamilyKey("Ascomycota", "MATsc")
SPANS = {"chr1": [1000, 5000]}
#: rec_c is a second MAT1-1 record, so withholding rec_a alone leaves a usable
#: reference and a miss is a genuine one.
FAMILIES = {"rec_a": MAT, "rec_b": MAT, "rec_c": MAT, "rec_sc": MATSC}
IDIOMORPHS = {"rec_a": {"MAT1-1"}, "rec_b": {"MAT1-2"}, "rec_c": {"MAT1-1"}, "rec_sc": {"MATa"}}


def _locus(contig="chr1", start=2000, end=3000, idiomorph="MAT1-1", family="Ascomycota:MAT"):
    return {"contig": contig, "start": start, "end": end, "idiomorph": idiomorph,
            "family": family, "confidence": "high"}


def _score(detected=(), suppressed=(), withheld=frozenset({"rec_a"}), record="rec_a"):
    return score_holdout(
        record_id=record, spans=SPANS, detected=list(detected), suppressed=list(suppressed),
        withheld=withheld, record_families=FAMILIES, record_idiomorphs=IDIOMORPHS,
    )


def test_right_place_right_idiomorph_is_a_hit():
    assert _score([_locus()]).status == "hit"


def test_right_place_wrong_idiomorph_is_not_a_hit():
    assert _score([_locus(idiomorph="MAT1-2")]).status == "wrong_idiomorph"


def test_an_undetermined_call_is_kept_apart():
    assert _score([_locus(idiomorph="undetermined")]).status == "hit_undetermined"


def test_a_correct_overlapping_locus_wins_over_a_wrong_one():
    s = _score([_locus(idiomorph="MAT1-2", start=1000, end=1500), _locus()])
    assert s.status == "hit"


def test_the_wrong_place_is_a_miss():
    s = _score([_locus(contig="chr2")])
    assert s.status == "miss"


def test_a_withheld_locus_at_the_right_place_is_a_bar_loss_not_a_miss():
    s = _score(suppressed=[_locus()])
    assert s.status == "suppressed"


def test_a_family_emptied_by_the_holdout_is_no_reference():
    s = _score(withheld=frozenset({"rec_a", "rec_b", "rec_c"}))
    assert s.status == "no_reference_family"


def test_only_the_other_idiomorph_left_is_no_reference():
    """Idiomorphs are non-homologous: MAT1-2 cannot find MAT1-1."""
    s = _score(withheld=frozenset({"rec_a", "rec_c"}))
    assert s.status == "no_reference_idiomorph"


def test_a_miss_with_a_usable_reference_is_a_genuine_miss():
    assert _score().status == "miss"


def test_a_combined_record_accepts_either_idiomorph():
    idio = dict(IDIOMORPHS, rec_a={"MAT1-1", "MAT1-2"})
    s = score_holdout(record_id="rec_a", spans=SPANS, detected=[_locus(idiomorph="MAT1-2")],
                      suppressed=[], withheld=frozenset({"rec_a"}),
                      record_families=FAMILIES, record_idiomorphs=idio)
    assert s.status == "hit"
