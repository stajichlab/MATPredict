"""Report-only cross-lineage idiomorph class (curator ruling 2026-10-01, B9).

Axis: MAT1-1 = alpha-box idiomorph, MAT1-2 = HMG idiomorph. Family names and
idiomorph labels do not change. S. pombe P and Yarrowia A/B are unassigned.
"""
from pathlib import Path

from MATPredict.detect.family_registry import FamilyKey, load_all_families
from MATPredict.detect.pipeline import idiomorph_class_for

DB = Path(__file__).resolve().parents[2] / "db"
FAMS = {f.key: f for f in load_all_families(DB)}


def cls(locus, label, phylum="Ascomycota"):
    return idiomorph_class_for(FAMS[FamilyKey(phylum, locus)], label)


def test_ascomycota_map():
    assert cls("MAT", "MAT1-1") == "MAT1-1" and cls("MATtub", "MAT1-2") == "MAT1-2"
    assert cls("MTL", "alpha") == "MAT1-1" and cls("MTL", "A") == "MAT1-2"
    assert cls("MATsc", "MATalpha") == "MAT1-1" and cls("MATsc", "MATa") == "MAT1-2"
    assert cls("PM", "M") == "MAT1-2" and cls("mat3", "M") == "MAT1-2"


def test_unsettled_labels_are_unassigned():
    assert cls("PM", "P") == "unassigned" and cls("mat2", "P") == "unassigned"
    assert cls("MATyl", "A") == "unassigned" and cls("MATyl", "B") == "unassigned"


def test_no_class_without_a_map_or_a_label():
    assert cls("MAT", "Plus", phylum="Mucoromycota") is None
    assert cls("HD", "A1", phylum="Basidiomycota") is None
    assert cls("MAT", "undetermined") is None and cls("MAT", None) is None


def test_every_ascomycota_family_maps_every_declared_idiomorph():
    for f in FAMS.values():
        if f.key.phylum == "Ascomycota":
            assert set(f.idiomorph_classes) == set(f.idiomorph_values or []), f.key
