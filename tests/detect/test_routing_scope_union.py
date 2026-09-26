"""A species-scoped family must not shadow a genus-scoped one.

`route` used to be a cascade that returned at the first tier matching anything:
exact-taxid matches first, lineage matches only if there were none. That let a
family scoped to a SPECIES suppress a family scoped to that species' GENUS.

Found by the leave-one-out benchmark on Schizosaccharomyces pombe (4896), as an
unexplained miss with a reference still present:

    mat2, mat3   taxonomic_scope [4896]   the SILENT cassettes
    mat1         taxonomic_scope [4895]   the ACTIVE locus (genus Schizosaccharomyces)

4895 is 4896's parent, so `mat1` should match by lineage -- but `direct` matched
mat2/mat3 and returned, so `mat1` was never searched. Detection looked like it
worked on S. pombe while never looking for the locus that determines mating type.

A genus-level scope is not weaker evidence than a species-level one; it is a
curator's statement about a different breadth.
"""
from MATPredict.detect.family_registry import Family, FamilyKey, route


def _fam(name, scope):
    return Family(key=FamilyKey("Ascomycota", name), vocabulary_type="enum",
                  idiomorph_values=["A", "B"], idiomorph_pattern=None,
                  genes=[{"name": "g1", "role": "core_MAT"}], taxonomic_scope=scope)


SPECIES, GENUS = 4896, 4895
ACTIVE = _fam("mat1", [GENUS])
SILENT_A, SILENT_B = _fam("mat2", [SPECIES]), _fam("mat3", [SPECIES])


def _lineage(taxid):
    return [SPECIES, GENUS, 4894] if taxid == SPECIES else [taxid]


def test_a_genus_scoped_family_survives_an_exact_species_match():
    d = route(SPECIES, [ACTIVE, SILENT_A, SILENT_B], lineage_taxids_resolver=_lineage)
    names = {f.key.locus_name for f in d.families}
    assert names == {"mat1", "mat2", "mat3"}, names


def test_an_exact_match_still_reports_direct():
    # The stronger signal stays visible; only the family SET changes.
    d = route(SPECIES, [ACTIVE, SILENT_A, SILENT_B], lineage_taxids_resolver=_lineage)
    assert d.routing_mode == "direct"


def test_lineage_only_is_still_lineage():
    d = route(SPECIES, [ACTIVE], lineage_taxids_resolver=_lineage)
    assert d.routing_mode == "lineage"
    assert {f.key.locus_name for f in d.families} == {"mat1"}


def test_no_family_is_listed_twice_when_both_tiers_match():
    both = _fam("both", [SPECIES, GENUS])
    d = route(SPECIES, [both], lineage_taxids_resolver=_lineage)
    assert [f.key.locus_name for f in d.families] == ["both"]


def test_an_unrelated_species_scoped_family_is_not_pulled_in():
    other = _fam("elsewhere", [999999])
    d = route(SPECIES, [ACTIVE, other], lineage_taxids_resolver=_lineage)
    assert {f.key.locus_name for f in d.families} == {"mat1"}
