"""A close vote between short fragments is reported `undetermined` (`thin_vote`, `--short-fragment-vote-policy`, #138)."""
import argparse
import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

from MATPredict.detect import cli as detect_cli
from MATPredict.detect import thin_vote as tv
from MATPredict.detect.pipeline import run_pipeline


def D(policy="undetermined", idiomorph="MAT1-1", margin=2.0, n=2, min_margin=5.0, covs=(40.0,)):
    return tv.decide(policy, idiomorph, margin, n, min_margin, list(covs))


def test_defaults_are_off():
    assert tv.DEFAULT_SHORT_FRAGMENT_VOTE_POLICY == "off"
    assert inspect.signature(run_pipeline).parameters["short_fragment_vote_policy"].default == "off"
    assert detect_cli._short_fragment_vote_policy_from_args(argparse.Namespace()) == "off"


def test_unknown_policy_is_an_error():
    with pytest.raises(ValueError):
        tv.check_policy("always")


def test_off_never_decides():
    assert D(policy="off") is None


def test_close_vote_on_short_fragments_is_undetermined_and_the_vote_is_kept():
    note = D(margin=0.7, covs=(35.0, 48.0))
    assert note["decision"] == "undetermined" and note["vote_idiomorph"] == "MAT1-1" and note["margin"] == 0.7
    assert note["best_core_coverage"] == 48.0 and note["min_margin"] == 5.0


@pytest.mark.parametrize("kw", [
    dict(covs=(75.0,)),                  # a fuller-length gene: a close vote is still a call
    dict(covs=(60.0,)),                  # exactly at the line stands
    dict(margin=5.0),                    # exactly at the minimum margin stands
    dict(margin=12.0),                   # a clear vote
    dict(n=1),                           # a single candidate has no vote to be close
    dict(margin=None),                   # no margin
    dict(idiomorph="undetermined"),      # already undetermined
    dict(idiomorph=None),
    dict(covs=()),                       # no coverage information: do not guess
    dict(covs=(None,)),
])
def test_everything_else_stands(kw):
    assert D(**kw) is None


def test_the_family_minimum_margin_is_the_line():
    assert D(margin=6.0, min_margin=8.0) is not None and D(margin=6.0, min_margin=5.0) is None


def test_cli_default_and_passthrough(monkeypatch, tmp_path):
    from MATPredict.detect.pipeline import DetectionOutcome

    monkeypatch.setattr(detect_cli, "MatpredictConfig", SimpleNamespace(from_env=lambda repo_root: SimpleNamespace(db_root=Path("db"))))
    seen: list = []

    def fake_run_pipeline(**kwargs):
        seen.append(kwargs.get("short_fragment_vote_policy", "MISSING"))
        return DetectionOutcome(results=[])

    monkeypatch.setattr(detect_cli, "build_reference_fasta", lambda db_root, out, family_keys=None, **kw: out)
    monkeypatch.setattr(detect_cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(detect_cli, "write_detection_gff3", lambda *a, **k: None)
    monkeypatch.setattr(detect_cli, "write_detection_report", lambda outcome, path, run=None: None)
    base = dict(genome="g.fa", proteins=None, taxid=None, out_dir=str(tmp_path / "out"), evidence_diagnostics=None, min_hits=1,
                min_identity=None, require_core_role=False, exclude_records="", emit_cds_fasta=False, phylum="Mucoromycota")
    assert detect_cli._cmd_detect(SimpleNamespace(**base)) == 0
    assert detect_cli._cmd_detect(SimpleNamespace(**base, short_fragment_vote_policy="undetermined")) == 0
    assert seen == ["off", "undetermined"]


def test_no_preset_sets_the_policy():
    # only the v0.9.1 preset turns this rule on (curator, 2026-11-07)
    assert {n for n, v in detect_cli.PRESETS.items() if "short_fragment_vote_policy" in v} == {"v0.9.1"}
    assert detect_cli.PRESETS["v0.9.1"]["short_fragment_vote_policy"] == "undetermined"
