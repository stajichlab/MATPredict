"""A second call the classifier could not assign, beside a determined call.

Curator's ruling 2026-09-27 (guard against spurious second calls). Adding the
tier-2 Umbelopsis records let weak flank hits (28-38%) be modelled next to HMG
paralogs, and 13 Mucoromycota genomes gained an 85-155 kb second call whose
idiomorph the HMM classifier scored but could not assign (margin 19.5-24.8,
below the 25-bit floor) while the genome already had a determined call of the
same family (results/2026-09-27_umbelopsis_curation/, results/2026-09-27_
mucoro_curation_guard/replay_guard.txt).

The rule is narrow on purpose. A replay showed a guard on every enum-vocabulary
family would withhold 8 Saccharomyces silent-cassette calls (MATALPHA2/MATA2
only), where several loci per genome are normal. So it applies only when the
classifier RAN on the call and returned `undetermined`.
"""
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.idiomorph import LOCUS_CLASS_MAT
from MATPredict.detect.pipeline import DetectionResult, GeneEvidence
from MATPredict.detect.secondary_undetermined import (
    WITHHELD_SECONDARY_UNDETERMINED, apply_secondary_undetermined_rule,
)

MUC = FamilyKey("Mucoromycota", "MAT")
OTHER = FamilyKey("Mucoromycota", "OTHER")


def _clf(idiomorph, margin):
    return {"method": "hmm", "idiomorph": idiomorph, "margin": margin,
            "scores": {"Minus": 50.0, "Plus": 50.0 - margin}, "min_margin": 25.0,
            "classifier_input": "model"}


def _call(contig, idiomorph, clf=None, key=MUC, locus_class=LOCUS_CLASS_MAT):
    ev = [GeneEvidence("sexM", "core_MAT", contig, 100, 400, "+", 40.0, 50.0,
                       "rec1", "tblastn_genome", status="polished_single")]
    return DetectionResult(
        family_key=key, contig=contig, start=100, end=400, confidence="medium",
        idiomorph=idiomorph, ambiguous_with=[], genes_found=["sexM"],
        genes_missing=[], fragmented=False, gene_evidence=ev,
        locus_class=locus_class, polished_genes=2, idiomorph_classifier=clf,
    )


def test_classifier_undetermined_beside_a_determined_call_is_withheld():
    first = _call("c1", "Minus", _clf("Minus", 97.6))
    second = _call("c2", "undetermined", _clf("undetermined", 24.4))
    kept, withheld = apply_secondary_undetermined_rule([first, second])
    assert kept == [first]
    assert [w.contig for w in withheld] == ["c2"]
    assert withheld[0].withheld_reason == WITHHELD_SECONDARY_UNDETERMINED


def test_a_lone_undetermined_call_is_kept():
    # Umbelopsis sp. AD052: its only call is undetermined -- nothing to defer to.
    only = _call("c1", "undetermined", _clf("undetermined", 23.9))
    kept, withheld = apply_secondary_undetermined_rule([only])
    assert kept == [only] and withheld == []


def test_undetermined_without_a_classifier_verdict_is_kept():
    # Saccharomyces cassettes: undetermined because only the non-voting
    # a2/alpha2 pair was found; no classifier ran.
    first = _call("c1", "a", None, key=OTHER)
    cassette = _call("c2", "undetermined", None, key=OTHER)
    kept, withheld = apply_secondary_undetermined_rule([first, cassette])
    assert kept == [first, cassette] and withheld == []


def test_homothallic_candidate_is_never_withheld():
    first = _call("c1", "Minus", _clf("Minus", 97.6))
    homo = _call("c2", "undetermined", _clf("undetermined", 5.0),
                 locus_class="homothallic_candidate")
    kept, withheld = apply_secondary_undetermined_rule([first, homo])
    assert homo in kept and withheld == []


def test_a_determined_call_in_another_family_does_not_count():
    other = _call("c1", "a", _clf("a", 80.0), key=OTHER)
    undet = _call("c2", "undetermined", _clf("undetermined", 20.0))
    kept, withheld = apply_secondary_undetermined_rule([other, undet])
    assert kept == [other, undet] and withheld == []


def test_two_determined_calls_are_untouched():
    a = _call("c1", "Minus", _clf("Minus", 97.6))
    b = _call("c2", "Plus", _clf("Plus", 60.0))
    kept, withheld = apply_secondary_undetermined_rule([a, b])
    assert kept == [a, b] and withheld == []
