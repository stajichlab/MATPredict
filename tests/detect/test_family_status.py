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


def test_cross_call_without_gene_evidence_has_equal_means_and_keeps_the_smaller_family_name():
    sc, mtl = _call(SC, 100, 5000, confidence="medium"), _call(MTL, 120, 5000, confidence="high")
    out = fs.resolve_cross_calls([sc, mtl], FAMILIES, {OTHER})
    assert [r.family_key for r in out] == [SC]
    assert out[0].cross_call[0]["family"] == "Ascomycota:MTL" and "near_tie" not in out[0].cross_call[0]
    assert out[0].family_ambiguous["kept_family"] == "Ascomycota:MATsc" and out[0].family_ambiguous["gap"] == 0.0


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


OG = FamilyKey("Ascomycota", "MATog")


def _og_families():
    f = _family(OG, [460523])
    from dataclasses import replace
    return FAMILIES + [replace(f, cassette_report=True)]


def _og(idiomorph, start, end, contig="c1"):
    from dataclasses import replace
    return replace(_call(OG, start, end, contig=contig), idiomorph=idiomorph)


def test_two_idiomorphs_are_indicated_as_cassettes_and_not_merged():
    a, alpha = _og("MATa", 1, 2000, "c1"), _og("MATalpha", 1, 3000, "c2")
    out = fs.apply([a, alpha], _og_families(), {460523}, "on")
    assert len(out) == 2
    assert out[0].cassette == {"idiomorph": "MATa", "partner_idiomorphs": ["MATalpha"], "role": "cassette_undetermined",
                               "note": out[0].cassette["note"]}
    assert out[1].cassette["partner_idiomorphs"] == ["MATa"]


def test_one_call_or_two_of_the_same_idiomorph_is_not_a_cassette():
    fams = _og_families()
    assert fs.apply([_og("MATa", 1, 2000)], fams, {460523}, "on")[0].cassette is None
    two = fs.apply([_og("MATa", 1, 2000, "c1"), _og("MATa", 1, 2000, "c2")], fams, {460523}, "on")
    assert all(r.cassette is None for r in two)


def test_a_family_without_cassette_report_is_never_indicated():
    out = fs.apply([_call(SC, 1, 100, contig="c1"), _call(SC, 1, 100, contig="c2")], FAMILIES, {OTHER}, "on")
    assert all(r.cassette is None for r in out)


def _scored(key, start, end, identity, confidence="high"):
    from dataclasses import replace
    c = _call(key, start, end, [_gene("G1", "core_MAT", start, end)], confidence=confidence)
    g = replace(c.gene_evidence[0], identity=identity)
    return replace(c, gene_evidence=[g])


def test_higher_core_identity_wins_not_the_tier():
    sc, mtl = _scored(SC, 100, 5000, 60.0, "medium"), _scored(MTL, 100, 5000, 50.0, "high")
    out = fs.resolve_cross_calls([sc, mtl], FAMILIES, {OTHER})
    assert [r.family_key for r in out] == [SC]
    assert out[0].cross_call[0]["best_core_identity"] == 50.0 and "near_tie" not in out[0].cross_call[0]


def test_there_is_no_margin_a_small_gap_is_decided_by_identity_and_flagged_ambiguous():
    sc, mtl = _scored(SC, 100, 5000, 48.0, "medium"), _scored(MTL, 100, 5000, 50.0, "high")
    out = fs.apply([sc, mtl], FAMILIES, {OTHER}, "on")
    assert [r.family_key for r in out] == [MTL]                  # 2 points higher wins; no preference for the family without a rule
    amb = out[0].family_ambiguous
    assert amb == {"reason": "cross_call_unscoped", "kept_family": "Ascomycota:MTL",
                   "identities": {"Ascomycota:MATsc": 48.0, "Ascomycota:MTL": 50.0}, "gap": 2.0,
                   "genome_mean_identity": {"Ascomycota:MATsc": 48.0, "Ascomycota:MTL": 50.0},
                   "genome_loci": {"Ascomycota:MATsc": 1, "Ascomycota:MTL": 1}}
    # the flank rule is a separate flag: MTL outside Serinales without PAP1, OBP1, PIK1 is unconfirmed
    assert out[0].family_status["status"] == fs.STATUS_UNCONFIRMED and out[0].family_status["reason"] == "mtl_outside_scope"


def test_exactly_equal_identity_keeps_the_smaller_family_name():
    sc, mtl = _scored(SC, 100, 5000, 50.0, "medium"), _scored(MTL, 100, 5000, 50.0, "high")
    assert [r.family_key for r in fs.resolve_cross_calls([sc, mtl], FAMILIES, {OTHER})] == [SC]


def test_the_family_is_decided_once_for_the_genome_by_the_mean_identity():
    # locus 1: MATsc 60 vs MTL 58 (per locus MATsc); locus 2: MATsc 40 vs MTL 55 (per locus MTL); means MATsc 50, MTL 56.5
    calls = [_scored(SC, 100, 5000, 60.0), _scored(MTL, 100, 5000, 58.0), _scored(SC, 100, 5000, 40.0), _scored(MTL, 100, 5000, 55.0)]
    from dataclasses import replace
    calls = [replace(c, contig="c1") if i < 2 else replace(c, contig="c2", gene_evidence=[replace(c.gene_evidence[0], contig="c2")]) for i, c in enumerate(calls)]
    out = fs.resolve_cross_calls(calls, FAMILIES, {OTHER})
    assert sorted((r.contig, r.family_key.locus_name) for r in out) == [("c1", "MTL"), ("c2", "MTL")]      # one family for the genome
    assert all(r.family_ambiguous["genome_mean_identity"] == {"Ascomycota:MATsc": 50.0, "Ascomycota:MTL": 56.5} for r in out)
    assert all(r.family_ambiguous["genome_loci"] == {"Ascomycota:MATsc": 2, "Ascomycota:MTL": 2} for r in out)


def test_a_scoped_decision_records_no_ambiguity_whatever_the_identities():
    from dataclasses import replace
    one = [_scored(SC, 100, 5000, 30.0), _scored(MTL, 100, 5000, 90.0)]
    two = [replace(c, contig="c2", gene_evidence=[replace(c.gene_evidence[0], contig="c2")]) for c in (_scored(SC, 100, 5000, 40.0), _scored(MTL, 100, 5000, 45.0))]
    out = fs.resolve_cross_calls(one + two, FAMILIES, {SACCHAROMYCETACEAE})
    assert sorted((r.contig, r.family_key.locus_name) for r in out) == [("c1", "MATsc"), ("c2", "MATsc")]
    assert all(r.family_ambiguous is None for r in out)


def test_a_third_family_group_is_decided_apart_from_the_scoped_group():
    from dataclasses import replace
    og = FamilyKey("Ascomycota", "MATog")
    fams = FAMILIES + [_family(og, [12345])]
    one = [_scored(SC, 100, 5000, 30.0), _scored(MTL, 100, 5000, 90.0)]                    # SC is scoped: MATsc wins, no ambiguity
    two = [replace(c, contig="c2", gene_evidence=[replace(c.gene_evidence[0], contig="c2")]) for c in (_scored(MTL, 100, 5000, 40.0), _scored(og, 100, 5000, 45.0))]
    out = {r.contig: r for r in fs.resolve_cross_calls(one + two, fams, {SACCHAROMYCETACEAE})}
    assert out["c1"].family_key == SC and out["c1"].family_ambiguous is None
    assert out["c2"].family_key == og and out["c2"].family_ambiguous["genome_loci"] == {"Ascomycota:MATog": 1, "Ascomycota:MTL": 1}   # c1 did not enter the means


def test_scope_still_beats_identity():
    sc, mtl = _scored(SC, 100, 5000, 30.0, "low"), _scored(MTL, 100, 5000, 90.0, "high")
    assert [r.family_key for r in fs.resolve_cross_calls([sc, mtl], FAMILIES, {SACCHAROMYCETACEAE})] == [SC]


# ---- the flank rule (curator, 2026-10-30): an unsupported out-of-scope MTL call cannot win against a family with its own flank gene
def _with_flank(call, name, role="flanking_conserved", start=None, end=None):
    from dataclasses import replace
    g = _gene(name, role, start if start is not None else call.start, end if end is not None else call.start + 100)
    return replace(call, gene_evidence=call.gene_evidence + [g], genes_found=call.genes_found + [name])


def test_an_mtl_call_without_flank_support_loses_to_a_family_with_its_own_flank_gene_whatever_the_identity():
    sc, mtl = _with_flank(_scored(SC, 100, 5000, 40.0, "medium"), "sla2"), _scored(MTL, 100, 5000, 70.0)
    out = fs.resolve_cross_calls([sc, mtl], FAMILIES, {OTHER})
    assert [r.family_key for r in out] == [SC]
    assert out[0].cross_call[0]["family"] == "Ascomycota:MTL" and out[0].cross_call[0]["flank_rule"] is True
    amb = out[0].family_ambiguous
    assert amb["reason"] == "cross_call_flank_rule" and amb["flank_rule_set_aside"] == ["Ascomycota:MTL"]
    assert amb["identities"] == {"Ascomycota:MATsc": 40.0, "Ascomycota:MTL": 70.0}


def test_an_mtl_call_with_two_flank_genes_keeps_the_identity_rule():
    mtl = _with_flank(_with_flank(_scored(MTL, 100, 5000, 70.0), "PAP1", "flanking_variable", 100, 200), "OBP1", "flanking_variable", 300, 400)
    sc = _with_flank(_scored(SC, 100, 5000, 40.0, "medium"), "sla2")
    out = fs.resolve_cross_calls([sc, mtl], FAMILIES, {OTHER})
    assert [r.family_key for r in out] == [MTL] and out[0].family_ambiguous["reason"] == "cross_call_unscoped"
    assert "flank_rule" not in out[0].cross_call[0]


def test_no_competitor_flank_gene_means_no_rule():
    sc, mtl = _scored(SC, 100, 5000, 40.0, "medium"), _scored(MTL, 100, 5000, 70.0)
    out = fs.resolve_cross_calls([sc, mtl], FAMILIES, {OTHER})
    assert [r.family_key for r in out] == [MTL] and out[0].family_ambiguous["reason"] == "cross_call_unscoped"


def test_an_mtl_call_inside_serinales_is_never_set_aside():
    sc, mtl = _with_flank(_scored(SC, 100, 5000, 40.0, "medium"), "sla2"), _scored(MTL, 100, 5000, 70.0)
    out = fs.resolve_cross_calls([sc, mtl], FAMILIES, {SERINALES})
    assert [r.family_key for r in out] == [MTL] and out[0].family_ambiguous is None


def test_a_set_aside_call_is_not_counted_in_the_genome_mean():
    from dataclasses import replace
    # locus 1: the rule sets MTL (90) aside; locus 2 is a plain identity decision. Means: MATsc (40+50)/2 = 45, MTL 55 over locus 2 only.
    l1 = [_with_flank(_scored(SC, 100, 5000, 40.0), "sla2"), _scored(MTL, 100, 5000, 90.0)]
    l2 = [replace(c, contig="c2", gene_evidence=[replace(g, contig="c2") for g in c.gene_evidence])
          for c in (_scored(SC, 100, 5000, 50.0), _scored(MTL, 100, 5000, 55.0))]
    out = fs.resolve_cross_calls(l1 + l2, FAMILIES, {OTHER})
    assert {r.contig: r.family_key.locus_name for r in out} == {"c1": "MATsc", "c2": "MTL"}
    c2 = next(r for r in out if r.contig == "c2")
    assert c2.family_ambiguous["genome_mean_identity"] == {"Ascomycota:MATsc": 45.0, "Ascomycota:MTL": 55.0}   # 90 is not in the MTL mean
