# tests/detect/test_idiomorph_policy.py
"""The optional-gene policy (`legacy` | `core_first`).

Parastagonospora pattern (v0.6.0, 87 loci): the HMG-like real gene MAT1-2-1
overlaps a MAT1-1-3 cross-hit, MAT1-1-3 scores 41.07% vs 40.0%, supersedes it,
and the locus is called MAT1-1. `core_first` keeps an accessory (`optional`)
gene from superseding or outvoting a core hit.
"""
from __future__ import annotations

import argparse

import pytest

from MATPredict.detect.cli import _idiomorph_policy_from_args
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.idiomorph import (
    DEFAULT_IDIOMORPH_POLICY, assign_idiomorph, idiomorph_candidates,
    resolve_idiomorph_overlaps,
)
from MATPredict.detect.search import SearchHit

FAM = Family(
    FamilyKey("Ascomycota", "MAT"), "enum", ["MAT1-1", "MAT1-2"], None,
    [
        {"name": "MAT1-1-1", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-1"]},
        {"name": "MAT1-2-1", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-2"]},
        {"name": "MAT1-1-3", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-1"], "optional": True},
        {"name": "MAT1-2-3", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-2"], "optional": True},
        {"name": "SLA2", "role": "flanking_conserved", "optional": True},
    ],
    [4890],
)


def _hit(gene, start, end, identity, bitscore=None, role="core_MAT"):
    return SearchHit(
        FAM.key, gene, role, "c1", start, end, "+", identity,
        f"rec_{gene}", "tblastn_genome", bitscore=bitscore,
    )


def _call(hits, policy):
    resolved, _ = resolve_idiomorph_overlaps(hits, FAM, policy=policy)
    live = [h for h in resolved if h.superseded_by is None]
    return assign_idiomorph(FAM, [h.gene_name for h in live], live, policy), resolved


def test_default_policy_is_legacy():
    assert DEFAULT_IDIOMORPH_POLICY == "legacy"
    assert _idiomorph_policy_from_args(argparse.Namespace()) == "legacy"
    hits = [_hit("MAT1-1-3", 100, 400, 41.07, 60.0), _hit("MAT1-2-1", 110, 410, 40.0, 58.0)]
    # omitting the argument is the same as legacy
    r_default, ev_default = resolve_idiomorph_overlaps(hits, FAM)
    r_legacy, ev_legacy = resolve_idiomorph_overlaps(hits, FAM, policy="legacy")
    assert r_default == r_legacy and ev_default == ev_legacy


def test_parastagonospora_pattern_legacy_picks_mat11_core_first_picks_mat12():
    hits = [_hit("MAT1-1-3", 100, 400, 41.07, 60.0), _hit("MAT1-2-1", 110, 410, 40.0, 58.0)]
    legacy, resolved = _call(hits, "legacy")
    assert legacy == "MAT1-1"
    assert {h.gene_name: h.superseded_by for h in resolved}["MAT1-2-1"] == "MAT1-1-3"
    core_first, resolved = _call(hits, "core_first")
    assert core_first == "MAT1-2"
    assert all(h.superseded_by is None for h in resolved)


def test_core_first_still_resolves_two_core_genes():
    # Botrytis-style disagreement between two core genes is NOT this bug and is
    # left to the existing first-pass rule under both policies.
    hits = [_hit("MAT1-1-1", 100, 400, 45.0, 70.0), _hit("MAT1-2-1", 110, 410, 40.0, 58.0)]
    for policy in ("legacy", "core_first"):
        call, resolved = _call(hits, policy)
        assert call == "MAT1-1"
        assert {h.gene_name: h.superseded_by for h in resolved}["MAT1-2-1"] == "MAT1-1-1"


def test_optional_only_locus_still_votes_under_core_first():
    hits = [_hit("MAT1-1-3", 100, 400, 50.0, 80.0), _hit("SLA2", 900, 1200, 70.0, 99.0, "flanking_conserved")]
    for policy in ("legacy", "core_first"):
        call, _ = _call(hits, policy)
        assert call == "MAT1-1"
    ranked = idiomorph_candidates(FAM, ["MAT1-1-3"], hits, "core_first")
    assert [c["idiomorph"] for c in ranked] == ["MAT1-1"]


def test_vote_ignores_optional_when_a_core_gene_has_a_hit():
    # Non-overlapping, so no resolution: only the vote decides.
    hits = [_hit("MAT1-2-3", 100, 300, 60.0, 120.0), _hit("MAT1-1-1", 5000, 5400, 40.0, 70.0)]
    assert assign_idiomorph(FAM, ["MAT1-2-3", "MAT1-1-1"], hits, "legacy") == "MAT1-2"
    assert assign_idiomorph(FAM, ["MAT1-2-3", "MAT1-1-1"], hits, "core_first") == "MAT1-1"


def test_unknown_policy_is_rejected():
    with pytest.raises(ValueError):
        resolve_idiomorph_overlaps([], FAM, policy="bogus")
