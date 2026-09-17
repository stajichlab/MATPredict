from __future__ import annotations
from pathlib import Path

import pytest

from MATPredict.detect.family_registry import Family, FamilyKey, load_all_families, route
from MATPredict.db.taxonomy import TaxonomyResult


def _family(phylum, name, scope):
    return Family(
        key=FamilyKey(phylum, name),
        vocabulary_type="enum",
        idiomorph_values=["a", "alpha"],
        idiomorph_pattern=None,
        genes=[{"name": "STE3", "role": "core_MAT"}],
        taxonomic_scope=scope,
    )


def test_route_narrows_to_matching_scope():
    families = [_family("Basidiomycota", "MAT", [5270]), _family("Ascomycota", "MATsc", [4930])]

    def fake_resolver(taxid):
        return TaxonomyResult(taxid=taxid, lineage="k__Fungi;...;g__Mycosarcoma;s__Mycosarcoma_maydis", is_current=True)

    result = route(5270, families, lineage_resolver=fake_resolver)
    assert [f.key.locus_name for f in result] == ["MAT"]


def test_route_falls_back_to_all_when_taxid_is_none():
    families = [_family("Basidiomycota", "MAT", [5270]), _family("Ascomycota", "MATsc", [4930])]
    assert route(None, families) == families


def test_route_falls_back_to_all_when_no_scope_matches():
    families = [_family("Basidiomycota", "MAT", [5270])]

    def fake_resolver(taxid):
        return TaxonomyResult(taxid=taxid, lineage="k__Fungi;...;g__Unrelated;s__Unrelated_sp", is_current=True)

    result = route(999999, families, lineage_resolver=fake_resolver)
    assert result == families


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
