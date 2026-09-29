"""The MAT-gene gate: is the classified core gene the MAT gene at all?

Curator's ruling 2026-09-28 (option a), on results/2026-09-28_validation_f3_f4/
NOTE.md: the classifier margin types sexM vs sexP reliably (0/540 wrong) but
does not separate MAT genes from HMG paralogs (46/189 paralog proteins pass
25 bits). An absolute score >= 100 bits on a MODELLED protein separates well
(96/108 true vs 9/189 paralogs); on fragments nothing does. So, for a family
with an idiomorph classifier:

* model-typed call whose best classifier score >= `mat_gene_min_score`
  (100): the core gene is the MAT gene -- kept;
* otherwise (a model score below it, or a fragment-typed call): kept only with
  flank support -- at least `flank_support_min_genes` (2) distinct roster
  flanking genes modelled at >= `flank_support_min_identity` (40%) in the call;
* otherwise withheld (`mat_gene_gate`).

Families without a classifier, calls without a classifier verdict, and
split-locus calls (which require flanks by construction) are untouched.
"""
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.mat_gene_gate import (
    DEFAULT_FLANK_SUPPORT_MIN_GENES, DEFAULT_FLANK_SUPPORT_MIN_IDENTITY,
    DEFAULT_MAT_GENE_MIN_SCORE, WITHHELD_MAT_GENE_GATE, apply_mat_gene_gate,
)
from MATPredict.detect.pipeline import DetectionResult, GeneEvidence

KEY = FamilyKey("Mucoromycota", "MAT")
OTHER = FamilyKey("Ascomycota", "MAT")
SPEC = {"type": "hmm", "dir": "/x", "min_margin": 25}


def _ev(gene, role, identity, status="polished_single", method="exonerate_refine"):
    return GeneEvidence(gene, role, "c1", 1, 100, "+", identity, None, "rec1",
                        method, status=status)


def _result(evidence, scores, classifier_input="model", idiomorph="Plus", key=KEY,
            split_locus=None):
    clf = None if scores is None else {
        "method": "hmm", "idiomorph": idiomorph, "scores": scores,
        "margin": abs(scores["Plus"] - scores["Minus"]),
        "classifier_input": classifier_input, "min_margin": 25,
    }
    return DetectionResult(
        family_key=key, contig="c1", start=1, end=1000, confidence="medium",
        idiomorph=idiomorph, ambiguous_with=[],
        genes_found=sorted({e.gene_name for e in evidence}), genes_missing=[],
        fragmented=False, gene_evidence=evidence, polished_genes=2,
        idiomorph_classifier=clf, split_locus=split_locus,
    )


CORE = _ev("sexP", "core_MAT", 38.0)
TWO_STRONG = [_ev("tptA", "flanking_conserved", 72.0), _ev("rnhA", "flanking_conserved", 84.0)]


def test_the_defaults_are_the_rulings_numbers():
    assert DEFAULT_MAT_GENE_MIN_SCORE == 100.0
    assert DEFAULT_FLANK_SUPPORT_MIN_IDENTITY == 40.0
    assert DEFAULT_FLANK_SUPPORT_MIN_GENES == 2


def test_a_model_scoring_100_bits_is_kept_without_flanks():
    r = _result([CORE], {"Plus": 157.0, "Minus": 35.6})
    kept, withheld = apply_mat_gene_gate([r], {KEY: SPEC})
    assert kept == [r] and withheld == []


def test_a_low_scoring_model_with_two_strong_flanks_is_kept():
    # Mucor irregularis shape: score 80.7, four flanks at 66-87%.
    r = _result([CORE] + TWO_STRONG, {"Plus": 50.6, "Minus": 80.7}, idiomorph="Minus")
    kept, withheld = apply_mat_gene_gate([r], {KEY: SPEC})
    assert kept == [r] and withheld == []


def test_a_low_scoring_model_with_one_strong_flank_is_withheld():
    # Syncephalastrum paralog shape: score 64.9, glrA 77% alone.
    r = _result([CORE, _ev("glrA", "flanking_variable", 77.2)],
                {"Plus": 55.3, "Minus": 64.9}, idiomorph="undetermined")
    kept, withheld = apply_mat_gene_gate([r], {KEY: SPEC})
    assert kept == []
    assert [w.withheld_reason for w in withheld] == [WITHHELD_MAT_GENE_GATE]
    assert withheld[0].withheld_detail["best_score"] == 64.9
    assert withheld[0].withheld_detail["supporting_flanks"] == ["glrA"]


def test_weak_or_unmodelled_flanks_do_not_count():
    # Circinella / Rhizomucor shape: flanks at 31-35%, one unmodelled.
    flanks = [_ev("glrA", "flanking_variable", 31.0),
              _ev("algA", "flanking_variable", 34.6),
              _ev("tptA", "flanking_conserved", 70.0, status="unpolished")]
    r = _result([CORE] + flanks, {"Plus": 36.8, "Minus": 17.7}, idiomorph="undetermined")
    kept, withheld = apply_mat_gene_gate([r], {KEY: SPEC})
    assert kept == [] and len(withheld) == 1


def test_the_same_flank_gene_twice_counts_once():
    flanks = [_ev("glrA", "flanking_variable", 77.0), _ev("glrA", "flanking_variable", 70.0)]
    r = _result([CORE] + flanks, {"Plus": 60.0, "Minus": 20.0})
    kept, withheld = apply_mat_gene_gate([r], {KEY: SPEC})
    assert kept == [] and len(withheld) == 1


def test_a_fragment_typed_call_needs_flank_support_even_at_high_score():
    # Absolute scores do not separate on fragments.
    r = _result([CORE], {"Plus": 121.8, "Minus": 29.6}, classifier_input="hsp_fragment")
    kept, withheld = apply_mat_gene_gate([r], {KEY: SPEC})
    assert kept == [] and len(withheld) == 1
    r2 = _result([CORE] + TWO_STRONG, {"Plus": 121.8, "Minus": 29.6},
                 classifier_input="hsp_fragment")
    kept, withheld = apply_mat_gene_gate([r2], {KEY: SPEC})
    assert kept == [r2] and withheld == []


def test_an_annotated_flank_counts_as_modelled():
    flanks = [_ev("tptA", "flanking_conserved", 72.0, status="not_polish_candidate",
                  method="diamond_proteome"),
              _ev("rnhA", "flanking_conserved", 84.0)]
    r = _result([CORE] + flanks, {"Plus": 60.0, "Minus": 20.0})
    kept, _ = apply_mat_gene_gate([r], {KEY: SPEC})
    assert kept == [r]


def test_families_without_a_classifier_or_verdict_are_untouched():
    r_other = _result([CORE], {"Plus": 30.0, "Minus": 20.0}, key=OTHER)
    r_none = _result([CORE], None)
    kept, withheld = apply_mat_gene_gate([r_other, r_none], {KEY: SPEC, OTHER: None})
    assert kept == [r_other, r_none] and withheld == []


def test_split_locus_calls_are_untouched():
    r = _result([CORE], {"Plus": 60.0, "Minus": 20.0}, split_locus={"core_gene": "sexP"})
    kept, withheld = apply_mat_gene_gate([r], {KEY: SPEC})
    assert kept == [r] and withheld == []


def test_the_thresholds_are_roster_overridable():
    spec = {**SPEC, "mat_gene_min_score": 50, "flank_support_min_genes": 1,
            "flank_support_min_identity": 70}
    r = _result([CORE], {"Plus": 60.0, "Minus": 20.0})
    assert apply_mat_gene_gate([r], {KEY: spec})[0] == [r]
    r2 = _result([CORE, _ev("glrA", "flanking_variable", 77.2)],
                 {"Plus": 45.0, "Minus": 20.0}, classifier_input="hsp_fragment")
    assert apply_mat_gene_gate([r2], {KEY: spec})[0] == [r2]
