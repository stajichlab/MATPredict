"""The S. pombe silent pair: a merged mat1 call holding a P and an M cluster about 11 kb apart is split into mat2 P and mat3 M (#254)."""
from __future__ import annotations

from dataclasses import replace

from MATPredict.detect import family_status as fs
from MATPredict.detect.family_registry import Family, FamilyKey, _gene_alias_map
from MATPredict.detect.pipeline import DetectionResult, GeneEvidence, LocusSegment
from MATPredict.detect.silent_pair import split_silent_pairs

MAT1, MAT2, MAT3 = (FamilyKey("Ascomycota", n) for n in ("mat1", "mat2", "mat3"))
POMBE, GENUS = 4896, 4895


def _fam(key, scope, genes, idio_class):
    return Family(key=key, vocabulary_type="enum", idiomorph_values=["M", "P"], idiomorph_pattern=None, genes=genes, taxonomic_scope=scope,
                  gene_aliases=_gene_alias_map(genes), idiomorph_classes=idio_class)


def _g(name, idio):
    return {"name": name, "role": "core_MAT", "present_in_idiomorphs": [idio]}


FAMILIES = [
    _fam(MAT1, [GENUS], [_g("mat1-Mc", "M"), _g("mat1-Mi", "M"), _g("mat1-Pc", "P"), _g("mat1-Pi", "P")], {"M": "MAT1-2", "P": "unassigned"}),
    _fam(MAT2, [POMBE], [_g("matPc", "P"), _g("matPi", "P")], {"P": "unassigned"}),
    _fam(MAT3, [POMBE], [_g("matMc", "M"), _g("matMi", "M")], {"M": "MAT1-2"}),
]


def _e(name, start, end, contig="c1"):
    return GeneEvidence(name, "core_MAT", contig, start, end, "+", 80.0, 100.0, "rec", "exonerate_refine", status="polished_agree")


def _merged(p_at, m_at, contig="c1"):
    ev = [_e("mat1-Pc", p_at, p_at + 353, contig), _e("mat1-Pi", p_at + 540, p_at + 1017, contig), _e("mat1-Mi", m_at, m_at + 125, contig), _e("mat1-Mc", m_at + 260, m_at + 802, contig)]
    lo, hi = min(e.start for e in ev), max(e.end for e in ev)
    return DetectionResult(family_key=MAT1, contig=contig, start=lo - 100, end=hi + 100, confidence="high", idiomorph="M", ambiguous_with=[],
                           genes_found=[e.gene_name for e in ev], genes_missing=[], fragmented=False, gene_evidence=ev, segments=[LocusSegment(contig, lo - 100, hi + 100)],
                           polished_genes=4)


def test_a_p_and_an_m_cluster_about_11_kb_apart_become_mat2_p_and_mat3_m():
    out = split_silent_pairs([_merged(329714, 342022)], FAMILIES, {POMBE})
    assert [(r.family_key.locus_name, r.idiomorph) for r in out] == [("mat2", "P"), ("mat3", "M")]
    p, m = out
    assert (p.start, p.end) == (329714, 330731) and (m.start, m.end) == (342022, 342824)
    assert p.genes_found == ["matPc", "matPi"] and m.genes_found == ["matMc", "matMi"]
    assert all(e.gene_name.startswith("mat") and not e.gene_name.startswith("mat1-") for e in p.gene_evidence + m.gene_evidence)
    assert p.cassette["role"] == "silent_pair" and p.cassette["partner_families"] == ["Ascomycota:mat3"] and m.cassette["partner_families"] == ["Ascomycota:mat2"]
    assert p.cassette["pair_distance_bp"] == 342022 - 330731 and p.confidence == m.confidence == "high" and p.polished_genes == m.polished_genes == 2
    assert p.idiomorph_class == "unassigned" and m.idiomorph_class == "MAT1-2" and p.genes_missing == m.genes_missing == []


def test_a_pair_farther_apart_than_20_kb_is_not_split_and_is_marked_ambiguous():
    out = split_silent_pairs([_merged(480787, 516366)], FAMILIES, {POMBE})
    assert len(out) == 1 and out[0].family_key == MAT1 and out[0].idiomorph == "M"
    assert out[0].cassette["role"] == "arrangement_ambiguous" and out[0].cassette["pair_distance_bp"] > 20_000


def test_a_lone_cassette_and_other_species_are_untouched():
    lone = DetectionResult(family_key=MAT1, contig="c9", start=1, end=1009, confidence="high", idiomorph="P", ambiguous_with=[], genes_found=["mat1-Pc", "mat1-Pi"],
                           genes_missing=[], fragmented=False, gene_evidence=[_e("mat1-Pc", 629, 982, "c9"), _e("mat1-Pi", 1, 477, "c9")])
    out = split_silent_pairs([_merged(329714, 342022), lone], FAMILIES, {POMBE})
    assert [(r.family_key.locus_name, r.idiomorph) for r in out] == [("mat2", "P"), ("mat3", "M"), ("mat1", "P")]      # the lone P stays mat1
    assert out[2].cassette is None
    # another Schizosaccharomyces species: mat2 and mat3 are not in its lineage, so nothing changes
    assert split_silent_pairs([_merged(329714, 342022)], FAMILIES, {GENUS, 4899})[0].family_key == MAT1


def test_the_family_status_rule_applies_it_and_never_changes_a_call_when_the_policy_is_off():
    merged = _merged(329714, 342022)
    assert [r.family_key.locus_name for r in fs.apply([merged], FAMILIES, {POMBE}, "on")] == ["mat2", "mat3"]
    assert fs.apply([merged], FAMILIES, {POMBE}, "off") == [merged]
