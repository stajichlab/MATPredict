"""A merged call is as cautious as its most cautious member (review F7).

results/2026-09-28_fable_review/: `replace(primary, ...)` kept only the
primary's verification, so a PR member's CAAX-unverified label was dropped;
and a low primary took the best member confidence even when members did not
carry the same idiomorph.
"""
from dataclasses import replace

from MATPredict.detect.locus_merge import merge_overlapping

from tests.detect.test_locus_merge import BA, GROUPS, PR, _r

UNV = {"status": "unverified", "reason": "CAAX-only in Agrocybaceae"}


def test_an_unverified_member_makes_the_merged_call_unverified():
    calls = [_r(BA, conf="high"), replace(_r(PR, conf="medium"), verification=UNV)]
    (m,) = merge_overlapping(calls, GROUPS)
    assert m.verification["status"] == "unverified"
    assert "CAAX-only in Agrocybaceae" in m.verification["reason"]


def test_reasons_of_several_unverified_members_are_combined():
    other = {"status": "unverified", "reason": "override route"}
    calls = [replace(_r(BA), verification=other), replace(_r(PR), verification=UNV)]
    (m,) = merge_overlapping(calls, GROUPS)
    assert "override route" in m.verification["reason"]
    assert "CAAX-only in Agrocybaceae" in m.verification["reason"]


def test_verified_members_leave_verification_empty():
    (m,) = merge_overlapping([_r(BA), _r(PR)], GROUPS)
    assert m.verification is None


def test_members_with_the_same_idiomorph_take_the_best_confidence():
    calls = [_r(BA, conf="low", idio="3"), _r(PR, conf="high", idio="3")]
    (m,) = merge_overlapping(calls, GROUPS)
    assert m.confidence == "high"


def test_a_low_primary_keeps_its_confidence_when_idiomorphs_differ():
    # compatible (one undetermined) but not the same: primary's confidence stands
    calls = [_r(BA, conf="low", idio="3"), _r(PR, conf="high", idio="undetermined")]
    (m,) = merge_overlapping(calls, GROUPS)
    assert m.family_key == BA
    assert m.confidence == "low"
