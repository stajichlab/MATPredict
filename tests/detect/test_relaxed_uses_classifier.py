"""A relaxed-pass call takes its idiomorph from the HMM classifier.

Found while measuring review finding F1 (2026-09-28): once the relaxed pass
could fire behind withheld strict candidates, 17 new Mucoromycota relaxed
calls came out labelled by the old bitscore vote with no classifier verdict,
because `_relaxed_results` never looked at the verdict that is computed for
every cluster.
"""
from MATPredict.detect.classifier import ClassifierVerdict
from MATPredict.detect.pipeline import EvidenceFloor, _relaxed_results

from tests.detect.test_relaxed_modelled_count import RICH, SEARCHABLE, _cluster


def _verdict(idiomorph, margin):
    return ClassifierVerdict(scores={"MAT1-1": 120.0, "MAT1-2": 120.0 - margin},
                             margin=margin, idiomorph=idiomorph, min_margin=25.0,
                             proteins_scored=1, genes_scored=["g1"])


def test_the_classifier_verdict_sets_the_relaxed_idiomorph():
    c = _cluster()
    (r,) = _relaxed_results([c], [RICH], searchable_genes=SEARCHABLE,
                            evidence_floor=EvidenceFloor(),
                            classifier_verdicts={(id(c), RICH.key): _verdict("MAT1-1", 60.0)})
    assert r.idiomorph == "MAT1-1"
    assert r.idiomorph_classifier["margin"] == 60.0
    assert r.idiomorph_candidates[0]["basis"] == "hmm_classifier"


def test_an_undetermined_verdict_leaves_the_relaxed_call_undetermined():
    c = _cluster()
    (r,) = _relaxed_results([c], [RICH], searchable_genes=SEARCHABLE,
                            evidence_floor=EvidenceFloor(),
                            classifier_verdicts={(id(c), RICH.key): _verdict("undetermined", 10.0)})
    assert r.idiomorph == "undetermined"


def test_without_a_verdict_the_old_rule_still_applies():
    c = _cluster()
    (r,) = _relaxed_results([c], [RICH], searchable_genes=SEARCHABLE, evidence_floor=EvidenceFloor())
    assert r.idiomorph_classifier is None
