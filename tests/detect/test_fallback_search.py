"""#228: a routed family that reports nothing because every hit was withheld also searches the phylum-fallback families."""
from dataclasses import replace

from MATPredict.detect import fallback_search as fb
from MATPredict.detect.family_registry import Family, FamilyKey, _gene_alias_map, fallback_families
from MATPredict.detect.pipeline import MODELLED_GENE_BAR, DetectionOutcome, DetectionResult

OG, SC, MTL, BASIDIO = (FamilyKey("Ascomycota", "MATog"), FamilyKey("Ascomycota", "MATsc"), FamilyKey("Ascomycota", "MTL"), FamilyKey("Basidiomycota", "MAT"))
OGATEA, SACCHAROMYCOTINA = 461281, 4891


def _family(key, scope):
    genes = [{"name": "X", "role": "core_MAT"}]
    return Family(key=key, vocabulary_type="enum", idiomorph_values=["a", "alpha"], idiomorph_pattern=None, genes=genes,
                  taxonomic_scope=scope, gene_aliases=_gene_alias_map(genes))


FAMILIES = [_family(OG, [OGATEA]), _family(SC, []), _family(MTL, []), _family(BASIDIO, [])]
LINEAGE = lambda taxid: [1, 4751, 4890, OGATEA]
PHYLUM = lambda taxid: "Ascomycota"


def _call(key, withheld=None):
    return DetectionResult(family_key=key, contig="c", start=1, end=100, confidence="high", idiomorph="alpha", ambiguous_with=[],
                           genes_found=["X"], genes_missing=[], fragmented=False, gene_evidence=[], withheld_reason=withheld)


def _outcome(results=(), suppressed=(), mode="lineage", attempted=(OG,)):
    return DetectionOutcome(results=list(results), suppressed_loci=list(suppressed), routing_mode=mode, families_attempted=list(attempted))


def test_the_fallback_families_are_the_phylum_ones_without_the_routed_family():
    d = fallback_families(53333, FAMILIES, [FAMILIES[0]], LINEAGE, PHYLUM)
    assert d.routing_mode == "phylum_fallback" and [f.key for f in d.families] == [SC, MTL]       # not MATog again, not the Basidiomycota family


def test_a_phylum_with_no_other_family_gives_nothing_to_add():
    d = fallback_families(53333, [FAMILIES[0], FAMILIES[3]], [FAMILIES[0]], LINEAGE, PHYLUM)
    assert d.families == [] and d.routing_mode == "not_searched"


def test_the_search_is_needed_only_when_nothing_was_reported_and_the_bar_withheld_a_hit():
    assert fb.needed(_outcome(suppressed=[_call(OG, MODELLED_GENE_BAR)]))
    assert not fb.needed(_outcome(results=[_call(OG)], suppressed=[_call(OG, MODELLED_GENE_BAR)]))      # a call stands
    assert not fb.needed(_outcome())                                                                    # nothing found at all
    assert not fb.needed(_outcome(suppressed=[_call(OG, "flank_carried")]))                              # withheld for another reason
    assert not fb.needed(_outcome(suppressed=[_call(OG, MODELLED_GENE_BAR)], mode="phylum_fallback"))   # the phylum was searched already


def test_merge_keeps_the_routed_mode_and_records_what_the_fallback_added():
    first = _outcome(suppressed=[_call(OG, MODELLED_GENE_BAR)])
    second = _outcome(results=[_call(MTL)], mode="phylum_fallback", attempted=(SC, MTL))
    out = fb.merge(first, second, "phylum_fallback")
    assert out.routing_mode == "lineage" and [r.family_key for r in out.results] == [MTL]
    assert out.families_attempted == [OG, SC, MTL] and len(out.suppressed_loci) == 1
    assert out.fallback_search["families"] == ["Ascomycota:MATsc", "Ascomycota:MTL"] and out.fallback_search["calls"] == 1
    assert replace(first).fallback_search is None                                                          # the first outcome is not changed
