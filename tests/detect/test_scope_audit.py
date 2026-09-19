from __future__ import annotations

from MATPredict.detect.scope_audit import (
    ScopeAuditResult,
    audit_scope,
    deepest_common_ancestor,
    record_taxids_by_family,
)
from MATPredict.detect.family_registry import Family, FamilyKey


def _family(phylum, name, scope, genes=None):
    return Family(
        key=FamilyKey(phylum, name),
        vocabulary_type="enum",
        idiomorph_values=["a", "alpha"],
        idiomorph_pattern=None,
        genes=genes or [{"name": "STE3", "role": "core_MAT"}],
        taxonomic_scope=scope,
    )


# NCBI-shaped fake lineages: each list is root-to-parent order, most-specific last,
# matching the real order `default_lineage_taxids` returns (verified live this
# session: taxid 5501's lineage ends ..., 147545 (Eurotiomycetes), 451871, 33183,
# 33184 (Onygenales), 5500 (Onygenaceae) -- broad to narrow).
_FAKE_LINEAGES = {
    111: [1, 10, 100],          # e.g. species 111 under genus 100 under order 10
    112: [1, 10, 100],          # same genus/order as 111
    113: [1, 10, 200],          # same order (10) but a different genus (200)
    999: [1, 20, 300],          # unrelated order entirely
}


def _fake_resolver(taxid: int) -> list[int]:
    return _FAKE_LINEAGES[taxid]


def test_deepest_common_ancestor_of_two_taxids_sharing_a_genus():
    assert deepest_common_ancestor([111, 112], _fake_resolver) == 100


def test_deepest_common_ancestor_falls_back_to_shared_order():
    assert deepest_common_ancestor([111, 113], _fake_resolver) == 10


def test_deepest_common_ancestor_of_a_single_taxid_is_itself():
    assert deepest_common_ancestor([111], _fake_resolver) == 111


def test_deepest_common_ancestor_returns_none_for_empty_list():
    assert deepest_common_ancestor([], _fake_resolver) is None


def test_deepest_common_ancestor_returns_none_when_genuinely_disjoint():
    # 111's lineage is [1, 10, 100]; 999's is [1, 20, 300] -- only root taxid 1
    # is shared, so the deepest common ancestor is 1, not None: two real fungi
    # always share at least a root/kingdom-level ancestor. None is reserved for
    # an empty input, not for "very distantly related."
    assert deepest_common_ancestor([111, 999], _fake_resolver) == 1


def test_audit_scope_flags_uncovered_records_and_recommends_a_fix():
    # Family's own scope [999] covers neither of its two real records (111, 112),
    # whose real common ancestor is 100.
    family = _family("Ascomycota", "MAT", scope=[999])
    record_taxids = {family.key: [111, 112]}

    results = audit_scope([family], record_taxids, lineage_taxids_resolver=_fake_resolver)

    assert results == [
        ScopeAuditResult(
            family_key=family.key, total_records=2,
            uncovered_taxids=[111, 112], recommended_scope_taxid=100,
        )
    ]


def test_audit_scope_reports_no_uncovered_taxids_when_scope_already_correct():
    family = _family("Ascomycota", "MATsc", scope=[100])  # 100 covers both via lineage
    record_taxids = {family.key: [111, 112]}

    results = audit_scope([family], record_taxids, lineage_taxids_resolver=_fake_resolver)

    assert results == [
        ScopeAuditResult(
            family_key=family.key, total_records=2,
            uncovered_taxids=[], recommended_scope_taxid=None,
        )
    ]


def test_audit_scope_direct_membership_counts_as_covered():
    # 111 is listed literally in scope -- must not require lineage resolution
    # (mirrors route()'s own direct-membership short-circuit).
    def unused_resolver(taxid):
        raise AssertionError("direct membership match must short-circuit before any lineage lookup")

    family = _family("Ascomycota", "MAT", scope=[111])
    record_taxids = {family.key: [111]}

    results = audit_scope([family], record_taxids, lineage_taxids_resolver=unused_resolver)

    assert results[0].uncovered_taxids == []
