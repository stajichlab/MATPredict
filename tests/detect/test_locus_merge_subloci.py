"""A merged A or B call carries its subloci as structured evidence.

Curator's ruling 2026-09-28, option (a), after the subloci literature review
(results/2026-09-28_subloci_literature/NOTE.md): Aalpha/Abeta (and
Balpha/Bbeta) are paralogous specificity units inside ONE mating-type locus,
so group A merges HD + Aalpha + Abeta into one A call where they overlap. The
sublocus detail is kept: each contributing family's label, genes, coordinates
and completeness.
"""
from dataclasses import replace

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.locus_merge import merge_overlapping

from tests.detect.test_locus_merge import BA, BB, HD, PR, _r

AA = FamilyKey("Basidiomycota", "Aalpha")
AB = FamilyKey("Basidiomycota", "Abeta")
GROUPS = {HD: ("A", True), AA: ("A", False), AB: ("A", False),
          PR: ("B", True), BA: ("B", False), BB: ("B", False)}


def test_hd_aalpha_abeta_overlapping_merge_into_one_a_call():
    (m,) = merge_overlapping([_r(HD, 1000, 20000), _r(AA, 1000, 9000), _r(AB, 8000, 20000)], GROUPS)
    assert {s["sublocus"] for s in m.subloci} == {"HD", "Aalpha", "Abeta"}


def test_each_sublocus_carries_genes_coordinates_and_completeness():
    aa = replace(_r(AA, 1000, 9000, genes=[("Z", 1000, 1500), ("Y", 2000, 2600)]), genes_missing=[])
    ab = replace(_r(AB, 8000, 20000, genes=[("Z", 9000, 9500)]), genes_missing=["Y"])
    (m,) = merge_overlapping([aa, ab, _r(HD, 1000, 20000)], GROUPS)
    by = {s["sublocus"]: s for s in m.subloci}
    assert by["Aalpha"]["genes"] == ["Z", "Y"]
    assert (by["Aalpha"]["contig"], by["Aalpha"]["start"], by["Aalpha"]["end"]) == ("c1", 1000, 9000)
    assert by["Aalpha"]["completeness"] == "complete"
    assert by["Abeta"]["completeness"] == "partial"
    assert by["Abeta"]["genes_missing"] == ["Y"]
    assert by["HD"]["generic"] is True and by["Aalpha"]["generic"] is False


def test_group_b_merge_carries_balpha_and_bbeta_subloci():
    (m,) = merge_overlapping([_r(PR), _r(BA), _r(BB)], GROUPS)
    assert {s["sublocus"] for s in m.subloci} == {"PR", "Balpha", "Bbeta"}


def test_an_unmerged_call_has_no_subloci():
    (r,) = merge_overlapping([_r(AA)], GROUPS)
    assert r.subloci == []


def test_the_roster_merges_aalpha_with_abeta():
    from pathlib import Path
    from MATPredict.detect.family_registry import load_all_families
    fams = {f.key.locus_name: f for f in load_all_families(Path(__file__).resolve().parents[2] / "db")
            if f.key.phylum == "Basidiomycota"}
    assert not fams["Aalpha"].merge_separately and not fams["Abeta"].merge_separately
