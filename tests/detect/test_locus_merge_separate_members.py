"""The `merge_separately` roster switch (default off).

The curator first asked (2026-09-28) for Aalpha and Abeta to stay separate,
then ruled after the subloci literature review that they merge (option a);
no roster locus sets the switch. These tests pin the switch's behaviour so it
can be turned on without code if the ruling changes.
"""
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.locus_merge import merge_overlapping

from tests.detect.test_locus_merge import BA, BB, HD, PR, _r

AA = FamilyKey("Basidiomycota", "Aalpha")
AB = FamilyKey("Basidiomycota", "Abeta")

SEP = {HD: ("A", True, False), AA: ("A", False, True), AB: ("A", False, True),
       PR: ("B", True, False), BA: ("B", False, False), BB: ("B", False, False)}
TOGETHER = {k: (g, gen, False) for k, (g, gen, _) in SEP.items()}


def _fams(results):
    return sorted(tuple(sorted(m["family"] for m in r.merged_from)) if r.merged_from
                  else (f"{r.family_key.phylum}:{r.family_key.locus_name}",) for r in results)


def test_aalpha_and_abeta_overlapping_stay_separate():
    out = merge_overlapping([_r(AA, 1000, 9000), _r(AB, 1000, 9000)], SEP)
    assert len(out) == 2


def test_hd_merges_with_one_sub_locus_not_both():
    out = merge_overlapping(
        [_r(HD, 1000, 20000), _r(AA, 1000, 9000), _r(AB, 12000, 20000)], SEP)
    assert len(out) == 2
    assert sum(1 for r in out if r.merged_from) == 1
    assert not any({m["family"] for m in r.merged_from} >= {"Basidiomycota:Aalpha", "Basidiomycota:Abeta"}
                   for r in out)


def test_with_the_setting_off_all_three_merge():
    out = merge_overlapping(
        [_r(HD, 1000, 20000), _r(AA, 1000, 9000), _r(AB, 1000, 20000)], TOGETHER)
    assert len(out) == 1


def test_group_b_is_unchanged():
    out = merge_overlapping([_r(PR), _r(BA), _r(BB)], SEP)
    assert len(out) == 1
    assert len(out[0].merged_from) == 3


def test_two_part_group_tuples_still_work():
    groups = {HD: ("A", True), AA: ("A", False)}
    assert len(merge_overlapping([_r(HD), _r(AA)], groups)) == 1

