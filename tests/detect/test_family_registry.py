from __future__ import annotations
from pathlib import Path

import pytest

from MATPredict.detect.family_registry import (
    DEFAULT_MAX_CLUSTER_GAP_BP,
    Family,
    FamilyKey,
    derive_max_cluster_gap,
    load_all_families,
    route,
)


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
    assert [f.key.locus_name for f in result.families] == ["MAT"]
    assert result.routing_mode == "direct"


def test_route_falls_back_to_all_when_taxid_is_none():
    families = [_family("Basidiomycota", "MAT", [5270]), _family("Ascomycota", "MATsc", [4930])]
    decision = route(None, families)
    assert decision.families == families
    assert decision.routing_mode == "exhaustive"


def test_route_falls_back_to_all_when_no_scope_matches():
    families = [_family("Basidiomycota", "MAT", [5270])]

    def fake_resolver(taxid):
        return [1, 2, 3]  # unrelated ancestor taxids, none of which is 5270

    result = route(
        999999, families,
        lineage_taxids_resolver=fake_resolver,
        phylum_name_resolver=lambda _taxid: None,  # phylum unresolvable -> no narrower fallback
    )
    assert result.families == families
    assert result.routing_mode == "exhaustive"


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
    assert [f.key.locus_name for f in result.families] == ["MAT"]
    assert result.routing_mode == "lineage"


def test_route_does_not_match_lineage_that_only_shares_unrelated_ancestors():
    families = [_family("Basidiomycota", "MAT", [5270]), _family("Ascomycota", "MATsc", [4930])]

    def fake_resolver(taxid):
        return [4751]  # shares kingdom Fungi with everything; not in either family's scope

    result = route(
        999999, families,
        lineage_taxids_resolver=fake_resolver,
        phylum_name_resolver=lambda _taxid: None,
    )
    # exhaustive fallback, not a false-positive kingdom-level match
    assert result.families == families
    assert result.routing_mode == "exhaustive"


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


def _three_phyla():
    """Two Ascomycota families and one Basidiomycota family, none of whose
    scopes will match the taxid used in the phylum-fallback tests below."""
    return [
        _family("Ascomycota", "MATsc", [4930]),
        _family("Ascomycota", "MATyl", [4951]),
        _family("Basidiomycota", "MAT", [5270]),
    ]


def test_route_falls_back_to_the_query_phylum_when_no_scope_matches():
    """The routing gap this fixes: taxid 294748's real lineage
    ([131567, 2759, 33154, 4751, 451864, 4890, 716545, 147537, 3239874,
    2916678, 766764, 5475, 5476]) intersects NO curated family's
    taxonomic_scope, so the old code returned all 19 families across all
    three phyla and the pipeline serially polished cross-phylum candidates.
    The lineage names 4890 (Ascomycota), so only Ascomycota's families are
    plausible."""
    families = _three_phyla()

    def fake_lineage(taxid):
        return [131567, 2759, 33154, 4751, 4890, 147537]

    result = route(
        294748, families,
        lineage_taxids_resolver=fake_lineage,
        phylum_name_resolver=lambda _taxid: "Ascomycota",
    )
    assert [f.key.locus_name for f in result.families] == ["MATsc", "MATyl"]
    assert result.routing_mode == "phylum_fallback"
    assert result.phylum == "Ascomycota"


def test_route_phylum_fallback_is_exhaustive_when_the_phylum_has_no_families():
    """A resolvable phylum with no curated families (e.g. Chytridiomycota,
    absent from db/) must not route to an EMPTY family set -- that would
    silently detect nothing. It degrades to the exhaustive set instead."""
    families = _three_phyla()
    result = route(
        999999, families,
        lineage_taxids_resolver=lambda _t: [4751],
        phylum_name_resolver=lambda _t: "Chytridiomycota",
    )
    assert result.families == families
    assert result.routing_mode == "exhaustive"


def test_route_phylum_resolver_failure_degrades_to_exhaustive():
    """A resolver failure (network error, unknown taxid) degrades the same
    way the lineage resolver's failure already does -- never raises out."""
    families = _three_phyla()

    def boom(_taxid):
        raise RuntimeError("NCBI unreachable")

    result = route(
        999999, families,
        lineage_taxids_resolver=lambda _t: [4751],
        phylum_name_resolver=boom,
    )
    assert result.families == families
    assert result.routing_mode == "exhaustive"


def test_route_explicit_phylum_skips_taxid_routing_entirely():
    """`matpredict detect --phylum Ascomycota` is an outright restriction:
    no direct check, no lineage fetch, no phylum fallback."""
    families = _three_phyla()

    def unused(_taxid):
        raise AssertionError("--phylum must skip taxid routing entirely")

    result = route(
        5270, families,
        lineage_taxids_resolver=unused,
        phylum_name_resolver=unused,
        phylum="Ascomycota",
    )
    assert [f.key.locus_name for f in result.families] == ["MATsc", "MATyl"]
    assert result.routing_mode == "explicit_phylum"
    assert result.phylum == "Ascomycota"


def test_route_explicit_phylum_applies_even_without_a_taxid():
    families = _three_phyla()
    result = route(None, families, phylum="Basidiomycota")
    assert [f.key.phylum for f in result.families] == ["Basidiomycota"]
    assert result.routing_mode == "explicit_phylum"


def test_available_phyla_is_read_from_db_root_at_runtime(tmp_path):
    """`--phylum`'s choices are discovered from db_root, never hardcoded, so
    a fourth phylum directory becomes selectable with no code change."""
    from MATPredict.detect.family_registry import available_phyla

    for phylum in ("Zoopagomycota", "Ascomycota"):
        (tmp_path / phylum).mkdir()
        (tmp_path / phylum / "order.yml").write_text(f"phylum: {phylum}\nloci: []\n")
    # Neither a candidates staging directory nor a schema directory declares an
    # order.yml, so neither can be offered as a --phylum choice.
    (tmp_path / "candidates").mkdir()
    (tmp_path / "_schema").mkdir()

    assert available_phyla(tmp_path) == ["Ascomycota", "Zoopagomycota"]


def test_available_phyla_agrees_with_the_families_actually_loaded(real_db_root):
    """A relationship, not a literal list: whatever phyla `db/` holds today,
    the `--phylum` choices must be exactly the phyla `load_all_families`
    produces `FamilyKey.phylum` values for. Asserting the three current names
    would break the moment curation adds a phylum, for a reason that has
    nothing to do with the routing code this guards."""
    from MATPredict.detect.family_registry import available_phyla

    phyla = available_phyla(real_db_root)
    assert phyla == sorted({f.key.phylum for f in load_all_families(real_db_root)})
    assert phyla, "the live database declares no phyla at all"


def test_available_phyla_skips_a_malformed_order_yml_and_warns(tmp_path, caplog):
    """A curator mid-edit can leave one `order.yml` unparseable. That must cost
    the CLI that phylum's `--phylum` choice and nothing else: `available_phyla`
    is called while the top-level argparse parser is being built, so raising
    here would abort `matpredict --help` and every unrelated subcommand too.
    The file is named in a warning so the breakage is not silent."""
    (tmp_path / "Ascomycota").mkdir()
    (tmp_path / "Ascomycota" / "order.yml").write_text("phylum: Ascomycota\nloci: []\n")
    (tmp_path / "Broken").mkdir()
    (tmp_path / "Broken" / "order.yml").write_text("phylum: [unclosed\n  - bad: :\n")

    from MATPredict.detect.family_registry import available_phyla

    with caplog.at_level("WARNING"):
        assert available_phyla(tmp_path) == ["Ascomycota"]
    assert "Broken/order.yml" in caplog.text


# --- per-locus cluster gap (order.yml `max_cluster_gap_bp`) ---


def _order_yml(tmp_path, phylum, extra_locus_lines=""):
    (tmp_path / phylum).mkdir(exist_ok=True)
    (tmp_path / phylum / "order.yml").write_text(
        f"phylum: {phylum}\nloci:\n  - locus_name: MAT\n    vocabulary_type: enum\n"
        "    idiomorph_values: [a, alpha]\n    taxonomic_scope: [1]\n"
        f"{extra_locus_lines}"
        "    genes:\n      - {name: STE3, role: core_MAT}\n"
    )


def test_load_all_families_reads_an_explicit_max_cluster_gap_bp(tmp_path):
    _order_yml(tmp_path, "Mucoro", "    max_cluster_gap_bp: 50000\n")
    families = load_all_families(tmp_path)
    assert [f.max_cluster_gap_bp for f in families] == [50_000]


def test_load_all_families_defaults_max_cluster_gap_bp_when_absent(tmp_path):
    _order_yml(tmp_path, "Asco")
    families = load_all_families(tmp_path)
    assert [f.max_cluster_gap_bp for f in families] == [DEFAULT_MAX_CLUSTER_GAP_BP]
    assert DEFAULT_MAX_CLUSTER_GAP_BP == 25_000


def test_derive_max_cluster_gap_takes_the_maximum_over_families():
    families = [
        _family("Asco", "MAT", [1]),
        Family(
            key=FamilyKey("Mucoro", "MAT"),
            vocabulary_type="enum",
            idiomorph_values=None,
            idiomorph_pattern=None,
            genes=[],
            taxonomic_scope=[2],
            max_cluster_gap_bp=50_000,
        ),
    ]
    assert derive_max_cluster_gap(families) == 50_000


def test_derive_max_cluster_gap_of_no_families_is_the_default():
    assert derive_max_cluster_gap([]) == DEFAULT_MAX_CLUSTER_GAP_BP
