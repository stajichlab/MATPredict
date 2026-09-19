"""Offline audit of every curated record's taxid against its own family's
declared taxonomic_scope -- the same direct-membership/lineage-overlap check
family_registry.route() uses at detection time, run against the whole curated
DB instead of one query taxid, so a stale scope value can be found and fixed
without running the detection pipeline at all."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import yaml

from MATPredict.db.taxonomy import default_lineage_taxids
from MATPredict.detect.family_registry import Family, FamilyKey


def record_taxids_by_family(db_root: Path) -> dict[FamilyKey, list[int]]:
    """Every accepted (non-candidate) record's own taxid, grouped by the
    FamilyKey it belongs to. Mirrors family_registry.load_record_families's
    glob and candidates/ exclusion, but collects taxids instead of an index."""
    by_family: dict[FamilyKey, list[int]] = {}
    for meta_path in sorted(db_root.glob("*/*/*/metadata.yaml")):
        if meta_path.relative_to(db_root).parts[0] == "candidates":
            continue
        doc = yaml.safe_load(meta_path.read_text())
        locus_name = (doc.get("mating_type") or {}).get("locus_name")
        taxid = (doc.get("taxonomy") or {}).get("taxid")
        if not locus_name or taxid is None:
            continue
        key = FamilyKey(meta_path.parents[2].name, locus_name)
        by_family.setdefault(key, []).append(taxid)
    return by_family


def deepest_common_ancestor(
    taxids: list[int], lineage_taxids_resolver: Callable[[int], list[int]]
) -> int | None:
    """The most specific taxid shared by every lineage in `taxids`, or the
    single taxid itself when only one is given, or None for an empty list.

    `lineage_taxids_resolver` (normally default_lineage_taxids) returns each
    taxid's ancestor chain broad-to-narrow (root first), never including the
    taxid itself -- so each taxid's own full lineage path, for this
    comparison, is `[*ancestors, taxid]`. Two real taxa always share at least
    a root-level ancestor, so this returns None only for an empty input, not
    for "very distantly related" -- callers that want to reject a too-broad
    recommendation should check the returned taxid's own rank/scientific name
    (e.g. via a live NCBI esummary lookup) before writing it into order.yml,
    the same live-verification discipline this project already requires
    everywhere else.
    """
    if not taxids:
        return None
    paths = [[*lineage_taxids_resolver(t), t] for t in taxids]
    # A taxonomic lineage is a genuine tree path (each taxid has exactly one
    # parent, never rejoining a sibling branch), so the deepest ancestor
    # shared by every path is exactly its longest common PREFIX -- position
    # by position, stop at the first position where the paths disagree.
    # zip(*paths) already stops at the shortest path's length.
    common: list[int] = []
    for level in zip(*paths):
        if len(set(level)) != 1:
            break
        common.append(level[0])
    return common[-1] if common else None


@dataclass(frozen=True)
class ScopeAuditResult:
    family_key: FamilyKey
    total_records: int
    uncovered_taxids: list[int]
    recommended_scope_taxid: int | None


def _covered(taxid: int, scope: list[int], ancestors: list[int]) -> bool:
    return taxid in scope or bool(set(ancestors) & set(scope))


def audit_scope(
    families: list[Family],
    record_taxids: dict[FamilyKey, list[int]],
    lineage_taxids_resolver: Callable[[int], list[int]] = default_lineage_taxids,
) -> list[ScopeAuditResult]:
    """For every family with at least one curated record, check whether
    route()'s own direct-membership-or-lineage-overlap test would cover each
    record's taxid, and recommend a single replacement scope taxid (the
    deepest common ancestor of every UNCOVERED record's taxid) when it would
    not. Families with zero uncovered records get `recommended_scope_taxid =
    None` -- there is nothing to fix.

    Deliberately does NOT wrap `lineage_taxids_resolver` calls in a
    try/except the way `family_registry.route()` does. `route()` degrades a
    resolver failure to a safe fallback (search everything); an audit tool
    has no such fallback, and swallowing a resolver failure here could
    recommend a bad scope value into order.yml on the strength of a
    transient network error rather than a real lineage. Raising loudly is
    the correct behavior here, not an oversight to be "fixed" into matching
    route()."""
    results = []
    for family in families:
        taxids = record_taxids.get(family.key, [])
        uncovered = []
        for taxid in taxids:
            ancestors = [] if taxid in family.taxonomic_scope else lineage_taxids_resolver(taxid)
            if not _covered(taxid, family.taxonomic_scope, ancestors):
                uncovered.append(taxid)
        recommended = (
            deepest_common_ancestor(uncovered, lineage_taxids_resolver) if uncovered else None
        )
        results.append(ScopeAuditResult(
            family_key=family.key, total_records=len(taxids),
            uncovered_taxids=uncovered, recommended_scope_taxid=recommended,
        ))
    return results
