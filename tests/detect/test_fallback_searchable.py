"""Scope-only families are never searched outside their scope.

Curator's ruling 2026-10-01: the regression check of curation-puccinio found
redPR "A1/A2 high" and wallMAT calls in phylum-fallback Agaricomycetes
(Rhizoctonia, Auricularia, Fomitiporia, Dacryopinax) -- STE3 homology, not
those lineages' loci. Families marked `fallback_searchable: false` are now
searched only for genomes inside their taxonomic scope.
"""
from pathlib import Path

from MATPredict.detect.family_registry import Family, FamilyKey, load_all_families, route

DB = Path(__file__).resolve().parents[2] / "db"


def fam(name, scope, fallback=True):
    return Family(key=FamilyKey("Basidiomycota", name), vocabulary_type="enum",
                  idiomorph_values=["A1"], idiomorph_pattern=None, genes=[],
                  taxonomic_scope=scope, fallback_searchable=fallback)


FAMS = [fam("PR", [5338]), fam("HD", [5338]), fam("redPR", [231213], fallback=False)]


def names(d):
    return sorted(f.key.locus_name for f in d.families)


def test_phylum_fallback_skips_scope_only_families():
    d = route(1, FAMS, lineage_taxids_resolver=lambda t: [1, 99],
              phylum_name_resolver=lambda t: "Basidiomycota")
    assert d.routing_mode == "phylum_fallback"
    assert names(d) == ["HD", "PR"]


def test_explicit_phylum_skips_scope_only_families():
    assert names(route(None, FAMS, phylum="Basidiomycota")) == ["HD", "PR"]


def test_in_scope_genome_still_searches_it():
    d = route(5, FAMS, lineage_taxids_resolver=lambda t: [5, 231213],
              phylum_name_resolver=lambda t: "Basidiomycota")
    assert d.routing_mode == "lineage" and names(d) == ["redPR"]


def test_exhaustive_still_searches_everything():
    d = route(None, FAMS, exhaustive=True)
    assert names(d) == ["HD", "PR", "redPR"]


def test_shipped_roster_flags_the_four_lineage_families():
    flagged = {f.key.locus_name for f in load_all_families(DB) if not f.fallback_searchable}
    assert flagged == {"redPR", "redHD", "rustHD", "wallMAT"}


def test_pr_scope_includes_boletales():
    pr = next(f for f in load_all_families(DB) if f.key == FamilyKey("Basidiomycota", "PR"))
    assert 68889 in pr.taxonomic_scope
