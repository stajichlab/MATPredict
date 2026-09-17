from __future__ import annotations

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.idiomorph import assign_idiomorph

ENUM_FAMILY = Family(
    FamilyKey("P", "MAT"), "enum", ["a", "alpha"], None,
    [
        {"name": "SXI1", "role": "core_MAT", "present_in_idiomorphs": ["alpha"]},
        {"name": "SXI2", "role": "core_MAT", "present_in_idiomorphs": ["a"]},
    ],
    [1],
)
PATTERN_FAMILY = Family(FamilyKey("P", "HD"), "pattern", None, "^A[0-9]+$",
                         [{"name": "HD1", "role": "core_MAT"}, {"name": "HD2", "role": "core_MAT"}], [1])


def test_enum_family_resolves_idiomorph_from_gene_presence():
    assert assign_idiomorph(ENUM_FAMILY, ["SXI1"]) == "alpha"
    assert assign_idiomorph(ENUM_FAMILY, ["SXI2"]) == "a"


def test_enum_family_undetermined_when_ambiguous_or_empty():
    assert assign_idiomorph(ENUM_FAMILY, []) == "undetermined"
    assert assign_idiomorph(ENUM_FAMILY, ["SXI1", "SXI2"]) == "undetermined"


def test_pattern_family_always_undetermined():
    assert assign_idiomorph(PATTERN_FAMILY, ["HD1", "HD2"]) == "undetermined"
