from __future__ import annotations
from pathlib import Path

import pytest

from MATPredict.detect.family_registry import Family, FamilyKey, load_all_families, route


def _family(phylum, name, scope):
    return Family(
        key=FamilyKey(phylum, name),
        vocabulary_type="enum",
        idiomorph_values=["a", "alpha"],
        idiomorph_pattern=None,
        genes=[{"name": "STE3", "role": "core_MAT"}],
        taxonomic_scope=scope,
    )


def test_route_narrows_to_matching_scope_by_direct_membership():
    families = [_family("Basidiomycota", "MAT", [5270]), _family("Ascomycota", "MATsc", [4930])]

    def unused_resolver(taxid):
        raise AssertionError("direct membership match must short-circuit before any lineage lookup")

    result = route(5270, families, lineage_taxids_resolver=unused_resolver)
    assert [f.key.locus_name for f in result] == ["MAT"]


def test_route_falls_back_to_all_when_taxid_is_none():
    families = [_family("Basidiomycota", "MAT", [5270]), _family("Ascomycota", "MATsc", [4930])]
    assert route(None, families) == families


def test_route_falls_back_to_all_when_no_scope_matches():
    families = [_family("Basidiomycota", "MAT", [5270])]

    def fake_resolver(taxid):
        return [1, 2, 3]  # unrelated ancestor taxids, none of which is 5270

    result = route(999999, families, lineage_taxids_resolver=fake_resolver)
    assert result == families


def test_route_matches_via_lineage_when_no_direct_membership():
    """A species taxid isn't itself in a family's broad taxonomic_scope, but its NCBI
    Taxonomy lineage includes the scope's subphylum-level taxid -- this is the 50/61
    curated-record routing gap the lineage-aware fix addresses."""
    # 5270 (Ustilago maydis, direct in this fixture) stands in for a subphylum/class
    # scope taxid (e.g. Ascomycota's real 222544 Pezizomycotina); 5271 stands in for a
    # descendant species/strain taxid that only appears in the family's scope via lineage.
    families = [_family("Basidiomycota", "MAT", [5270]), _family("Ascomycota", "MATsc", [4930])]

    def fake_resolver(taxid):
        assert taxid == 5271
        return [4751, 5270]  # e.g. k__Fungi (4751) -> ... -> the broad scope taxid (5270)

    result = route(5271, families, lineage_taxids_resolver=fake_resolver)
    assert [f.key.locus_name for f in result] == ["MAT"]


def test_route_does_not_match_lineage_that_only_shares_unrelated_ancestors():
    families = [_family("Basidiomycota", "MAT", [5270]), _family("Ascomycota", "MATsc", [4930])]

    def fake_resolver(taxid):
        return [4751]  # shares kingdom Fungi with everything; not in either family's scope

    result = route(999999, families, lineage_taxids_resolver=fake_resolver)
    assert result == families  # exhaustive fallback, not a false-positive kingdom-level match


def test_load_all_families_reads_real_order_yml(tmp_path):
    (tmp_path / "TestPhylum").mkdir()
    (tmp_path / "TestPhylum" / "order.yml").write_text(
        "phylum: TestPhylum\n"
        "loci:\n"
        "  - locus_name: MAT\n"
        "    vocabulary_type: enum\n"
        "    idiomorph_values: [a, alpha]\n"
        "    taxonomic_scope: [4930]\n"
        "    genes:\n"
        "      - {name: STE3, role: core_MAT}\n"
    )
    families = load_all_families(tmp_path)
    assert families == [_family("TestPhylum", "MAT", [4930])]
