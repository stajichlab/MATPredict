"""The single-candidate identity gap reported apart from the vote margin (`margin_basis`, `--idiomorph-margin-policy`, #145)."""
import argparse
import inspect
from dataclasses import replace
from types import SimpleNamespace

import pytest

from MATPredict.detect import cli as detect_cli
from MATPredict.detect import margin_basis as mb
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import DetectionResult, run_pipeline

KEY = FamilyKey("Ascomycota", "MAT")
RES = [{"winner": "MAT1-2-1", "loser": "MAT1-1-3", "margin": 4.3}]


def _call(margin, candidates, resolutions=None, classifier=None):
    return DetectionResult(
        family_key=KEY, contig="c1", start=1, end=100, confidence="high", idiomorph="MAT1-2", ambiguous_with=[],
        genes_found=["MAT1-2-1"], genes_missing=[], fragmented=False, gene_evidence=[],
        idiomorph_candidates=candidates, idiomorph_margin=margin, idiomorph_resolutions=resolutions or [], idiomorph_classifier=classifier)


def test_defaults_are_off_everywhere():
    assert mb.DEFAULT_IDIOMORPH_MARGIN_POLICY == "off"
    assert inspect.signature(run_pipeline).parameters["idiomorph_margin_policy"].default == "off"
    assert detect_cli._idiomorph_margin_policy_from_args(argparse.Namespace()) == "off"


def test_unknown_policy_is_an_error():
    with pytest.raises(ValueError):
        mb.check_policy("on")


def test_off_returns_the_results_untouched():
    r = _call(4.3, [{"idiomorph": "MAT1-2", "score": 195.0}], RES)
    assert mb.apply([r], "off") == [r]


def test_separate_moves_a_single_candidate_fallback_margin():
    r = _call(4.3, [{"idiomorph": "MAT1-2", "score": 195.0}], RES)
    (out,) = mb.apply([r], "separate")
    assert out.idiomorph_margin is None and out.idiomorph_identity_gap == 4.3
    assert (out.idiomorph, out.confidence, out.start, out.end, out.locus_class) == (r.idiomorph, r.confidence, r.start, r.end, r.locus_class)


@pytest.mark.parametrize("r", [
    _call(12.0, [{"idiomorph": "MAT1-1", "score": 90.0}, {"idiomorph": "MAT1-2", "score": 78.0}], RES),      # a real vote
    _call(6.5, [{"idiomorph": "MAT1-2", "score": 90.0}], RES, classifier={"scores": {"MAT1-2": 9.0}}),       # the HMM classifier's margin
    _call(None, [{"idiomorph": "MAT1-2", "score": 90.0}], RES),                                              # no margin at all
    _call(4.3, [{"idiomorph": "MAT1-2", "score": 90.0}], []),                                                # no resolution behind it
])
def test_separate_never_touches_a_vote_classifier_or_empty_margin(r):
    assert mb.apply([r], "separate") == [r]


def test_the_report_writes_the_gap_only_when_set():
    from MATPredict.detect import report
    base = _call(4.3, [{"idiomorph": "MAT1-2", "score": 195.0}], RES)
    (moved,) = mb.apply([base], "separate")
    src = inspect.getsource(report)
    assert '"idiomorph_identity_gap"' in src
    assert moved.idiomorph_identity_gap == 4.3 and base.idiomorph_identity_gap is None


def test_cli_default_and_passthrough(monkeypatch, tmp_path):
    from pathlib import Path
    from MATPredict.detect.pipeline import DetectionOutcome

    monkeypatch.setattr(detect_cli, "MatpredictConfig", SimpleNamespace(from_env=lambda repo_root: SimpleNamespace(db_root=Path("db"))))
    seen: list = []

    def fake_run_pipeline(**kwargs):
        seen.append(kwargs.get("idiomorph_margin_policy", "MISSING"))
        return DetectionOutcome(results=[])

    monkeypatch.setattr(detect_cli, "build_reference_fasta", lambda db_root, out, family_keys=None, **kw: out)
    monkeypatch.setattr(detect_cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(detect_cli, "write_detection_gff3", lambda *a, **k: None)
    monkeypatch.setattr(detect_cli, "write_detection_report", lambda outcome, path, run=None: None)
    base = dict(genome="g.fa", proteins=None, taxid=None, out_dir=str(tmp_path / "out"), evidence_diagnostics=None, min_hits=1,
                min_identity=None, require_core_role=False, exclude_records="", emit_cds_fasta=False, phylum="Mucoromycota")
    assert detect_cli._cmd_detect(SimpleNamespace(**base)) == 0
    assert detect_cli._cmd_detect(SimpleNamespace(**base, idiomorph_margin_policy="separate")) == 0
    assert seen == ["off", "separate"]


def test_no_preset_sets_the_margin_policy():
    assert all("idiomorph_margin_policy" not in v for v in detect_cli.PRESETS.values())
