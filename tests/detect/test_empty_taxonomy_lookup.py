# tests/detect/test_empty_taxonomy_lookup.py
"""An EMPTY taxonomy answer is recorded in `routing_error`, not left blank (#324).

Seven *C. fimbriata* assemblies (taxid 5158, a valid taxon) were `not_searched` in the v0.9.0 run with
`routing_error: None`: the lineage and phylum lookups returned nothing (the shared cache had no document for 5158 and the
compute node had no network), which raises no exception, so nothing was recorded. The routing is unchanged (still
`not_searched`); only the report field and the not-searched reason now say that the lookup was empty.
"""
from MATPredict.detect.family_registry import Family, FamilyKey, route


def _fam(name, scope):
    return Family(key=FamilyKey("Ascomycota", name), vocabulary_type="enum",
                  idiomorph_values=["A", "B"], idiomorph_pattern=None,
                  genes=[{"name": "g1", "role": "core_MAT"}], taxonomic_scope=scope)


FAMS = [_fam("MAT", [4890])]


def test_an_empty_answer_from_both_lookups_is_recorded():
    d = route(5158, FAMS, lineage_taxids_resolver=lambda t: [], phylum_name_resolver=lambda t: None)
    assert d.routing_mode == "not_searched"
    assert d.routing_error and "5158" in d.routing_error and "no lineage and no phylum" in d.routing_error


def test_a_lineage_without_a_phylum_rank_is_not_an_error():
    d = route(5158, FAMS, lineage_taxids_resolver=lambda t: [131567, 2759], phylum_name_resolver=lambda t: None)
    assert d.routing_mode == "not_searched"
    assert d.routing_error is None


def test_an_exception_keeps_its_own_message():
    def boom(t):
        raise OSError("no network")
    d = route(5158, FAMS, lineage_taxids_resolver=boom, phylum_name_resolver=boom)
    assert "lineage lookup failed" in d.routing_error and "no lineage and no phylum" not in d.routing_error


def test_a_resolved_taxid_is_unchanged():
    d = route(5158, FAMS, lineage_taxids_resolver=lambda t: [4890], phylum_name_resolver=lambda t: "Ascomycota")
    assert d.routing_mode == "lineage" and d.routing_error is None
