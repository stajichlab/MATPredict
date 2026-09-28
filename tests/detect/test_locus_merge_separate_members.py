"""Aalpha and Abeta never merge with each other (curator, 2026-09-28, PROVISIONAL).

Group A keeps HD (generic) + Aalpha + Abeta, but each sub-locus merges only
with the generic HD call; Aalpha and Abeta stay separate calls, even when one
HD call overlaps both. Set per roster locus (`merge_separately: true`), so it
can be flipped without code while the subloci literature review runs
(results/2026-09-28_subloci_literature/). Group B (PR + Balpha + Bbeta) does
not set it and still merges into one call.
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


def test_the_roster_sets_it_for_aalpha_and_abeta_only(tmp_path):
    from pathlib import Path
    from MATPredict.detect.family_registry import load_all_families as load_families
    fams = {f.key.locus_name: f for f in load_families(Path(__file__).resolve().parents[2] / "db")
            if f.key.phylum == "Basidiomycota"}
    assert fams["Aalpha"].merge_separately and fams["Abeta"].merge_separately
    assert not any(fams[n].merge_separately for n in ("HD", "PR", "Balpha", "Bbeta") if n in fams)
