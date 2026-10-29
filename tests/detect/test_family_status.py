import pytest

from MATPredict.detect import family_status as fs
from MATPredict.detect.family_registry import Family, FamilyKey, _gene_alias_map
from MATPredict.detect.pipeline import DetectionResult, GeneEvidence

SC, MTL = FamilyKey("Ascomycota", "MATsc"), FamilyKey("Ascomycota", "MTL")
SERINALES, SACCHAROMYCETACEAE, OTHER = 2916678, 4893, 999


def _family(key, scope, rule=None):
    genes = [{"name": "X", "role": "core_MAT"}]
    return Family(
        key=key, vocabulary_type="enum", idiomorph_values=["a", "alpha"], idiomorph_pattern=None, genes=genes,
        taxonomic_scope=scope, gene_aliases=_gene_alias_map(genes), confirm_outside_scope=rule)


MTL_RULE = {"flank_genes": ["PAP1", "OBP1", "PIK1"], "min_flank_genes": 2, "window_bp": 15000}
FAMILIES = [_family(SC, [SACCHAROMYCETACEAE]), _family(MTL, [SERINALES], MTL_RULE)]


def _gene(name, role, start, end, status="polished_agree", contig="c1"):
    return GeneEvidence(gene_name=name, role=role, contig=contig, start=start, end=end, strand="+", identity=50.0,
                        coverage=50.0, reference_record_id="r", method="m", status=status)


def _call(key, start, end, genes=(), confidence="high", contig="c1", genes_found=None):
    found = genes_found if genes_found is not None else [g.gene_name for g in genes]
    return DetectionResult(family_key=key, contig=contig, start=start, end=end, confidence=confidence, idiomorph="A",
                           ambiguous_with=[], genes_found=found, genes_missing=[], fragmented=False,
                           gene_evidence=list(genes))


def test_off_returns_the_calls_untouched():
    calls = [_call(SC, 100, 5000), _call(MTL, 100, 5000)]
    assert fs.apply(calls, FAMILIES, {OTHER}, "off") == calls


def test_unknown_policy_is_an_error():
    with pytest.raises(ValueError):
        fs.check_policy("maybe")


def test_cross_call_keeps_the_family_scoped_to_the_genome():
    sc, mtl = _call(SC, 100, 5000, confidence="low"), _call(MTL, 100, 5000, confidence="high")
    out = fs.resolve_cross_calls([sc, mtl], FAMILIES, {SACCHAROMYCETACEAE})
    assert [r.family_key for r in out] == [SC]          # scoped wins even though MTL scores better
    assert out[0].cross_call[0]["family"] == "Ascomycota:MTL"


def test_cross_call_neither_scoped_keeps_the_better_score():
    sc, mtl = _call(SC, 100, 5000, confidence="medium"), _call(MTL, 120, 5000, confidence="high")
    out = fs.resolve_cross_calls([sc, mtl], FAMILIES, {OTHER})
    assert [r.family_key for r in out] == [MTL]
    assert out[0].cross_call[0]["family"] == "Ascomycota:MATsc"


def test_cross_call_needs_90_percent_overlap_on_the_same_contig():
    a, b = _call(SC, 100, 1000), _call(MTL, 900, 5000)
    assert len(fs.resolve_cross_calls([a, b], FAMILIES, {OTHER})) == 2
    c = _call(MTL, 100, 1000, contig="c2")
    assert len(fs.resolve_cross_calls([a, c], FAMILIES, {OTHER})) == 2


def test_same_family_overlaps_are_not_cross_calls():
    a, b = _call(SC, 100, 1000), _call(SC, 100, 1000)
    assert len(fs.resolve_cross_calls([a, b], FAMILIES, {OTHER})) == 2


def _mtl_with_flanks(names):
    core = [_gene("MTLA1", "core_MAT", 5000, 5200), _gene("MTLA2", "core_MAT", 4000, 4300)]
    flank = {"PAP1": (1000, 2500), "OBP1": (2700, 3900), "PIK1": (30000, 33000)}  # PIK1 is 25 kb away
    return _call(MTL, 1000, 5200, core + [_gene(n, "flanking_variable", *flank[n]) for n in names])


def test_outside_scope_with_two_flanks_in_the_window_is_confirmed():
    out = fs.flag_family_status([_mtl_with_flanks(["PAP1", "OBP1"])], FAMILIES, {OTHER})
    st = out[0].family_status
    assert st["status"] == fs.STATUS_CONFIRMED and st["reason"] == "mtl_outside_scope"
    assert st["flank_genes_found"] == ["OBP1", "PAP1"]


def test_a_flank_beyond_the_window_does_not_count():
    out = fs.flag_family_status([_mtl_with_flanks(["PAP1", "PIK1"])], FAMILIES, {OTHER})
    assert out[0].family_status["status"] == fs.STATUS_UNCONFIRMED
    assert out[0].family_status["flank_genes_found"] == ["PAP1"]


def test_no_flank_is_unconfirmed_and_confidence_is_unchanged():
    out = fs.flag_family_status([_mtl_with_flanks([])], FAMILIES, {OTHER})
    assert out[0].family_status["status"] == fs.STATUS_UNCONFIRMED
    assert out[0].confidence == "high"


def test_inside_scope_or_unknown_lineage_is_never_flagged():
    call = _mtl_with_flanks([])
    assert fs.flag_family_status([call], FAMILIES, {SERINALES})[0].family_status is None
    assert fs.flag_family_status([call], FAMILIES, set())[0].family_status is None


def test_family_without_a_rule_is_never_flagged():
    assert fs.flag_family_status([_call(SC, 1, 100)], FAMILIES, {OTHER})[0].family_status is None


def test_apply_resolves_then_flags():
    sc = _call(SC, 1000, 5200, confidence="low")
    mtl = _mtl_with_flanks(["PAP1", "OBP1"])
    out = fs.apply([sc, mtl], FAMILIES, {OTHER}, "on")
    assert [r.family_key for r in out] == [MTL]
    assert out[0].family_status["status"] == fs.STATUS_CONFIRMED
    assert out[0].cross_call[0]["family"] == "Ascomycota:MATsc"
