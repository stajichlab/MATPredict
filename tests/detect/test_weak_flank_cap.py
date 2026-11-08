"""The weak-flank tier cap (`weak_flank_cap`, `--weak-flank-cap-policy`, #129)."""
import argparse
import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

from MATPredict.detect import cli as detect_cli
from MATPredict.detect import weak_flank_cap as wf
from MATPredict.detect.pipeline import run_pipeline

FAM = ("Ascomycota", "MAT")


def S(gene, role, status="unpolished", identity=40.0, coverage=10.0):
    return wf.FlankStat(gene, role, status, identity, coverage)


CORE = S("MAT1-1-1", "core_MAT", "polished_agree", 80.0, 95.0)


def test_defaults_are_off():
    assert wf.DEFAULT_WEAK_FLANK_CAP_POLICY == "off"
    assert inspect.signature(run_pipeline).parameters["weak_flank_cap_policy"].default == "off"
    assert detect_cli._weak_flank_cap_policy_from_args(argparse.Namespace()) == "off"


def test_unknown_policy_is_an_error():
    with pytest.raises(ValueError):
        wf.check_policy("maybe")


def test_off_never_caps():
    assert wf.evaluate("high", "MAT1-1", FAM, [CORE, S("APN2", "flanking_conserved")], policy="off") is None


def test_all_weak_flank_genes_cap_a_high_call():
    note = wf.evaluate("high", "MAT1-1", FAM, [CORE, S("APN2", "flanking_conserved", coverage=5.0)], policy="on")
    assert note["from"] == "high" and note["to"] == "medium" and [g["gene"] for g in note["weak_flank_genes"]] == ["APN2"]


@pytest.mark.parametrize("flank", [
    S("APN2", "flanking_conserved", status="polished_agree", identity=30.0, coverage=10.0),   # polished model
    S("APN2", "flanking_conserved", identity=52.0, coverage=31.0),                            # unpolished but strong enough
])
def test_one_supported_flank_gene_keeps_the_tier(flank):
    assert wf.evaluate("high", "MAT1-1", FAM, [CORE, S("SLA2", "flanking_conserved"), flank], policy="on") is None


def test_identity_alone_or_coverage_alone_is_not_support():
    assert wf.evaluate("high", "MAT1-1", FAM, [CORE, S("APN2", "flanking_conserved", identity=70.0, coverage=10.0)], policy="on")
    assert wf.evaluate("high", "MAT1-1", FAM, [CORE, S("APN2", "flanking_conserved", identity=30.0, coverage=80.0)], policy="on")


@pytest.mark.parametrize("confidence,idiomorph,family,stats", [
    ("medium", "MAT1-1", FAM, [CORE, S("APN2", "flanking_conserved")]),             # only high calls are touched
    ("high", "MAT1-1", FAM, [CORE]),                                                # no flank gene: nothing weak lifted it
    ("high", "undetermined", FAM, [CORE, S("APN2", "flanking_conserved")]),         # idiomorph not called
    ("high", "MAT1-1", ("Ascomycota", "MATsc"), [CORE, S("SLA2", "flanking_conserved")]),  # outside the measured family
])
def test_out_of_scope_calls_are_untouched(confidence, idiomorph, family, stats):
    assert wf.evaluate(confidence, idiomorph, family, stats, policy="on") is None


def test_cli_default_and_passthrough(monkeypatch, tmp_path):
    from MATPredict.detect.pipeline import DetectionOutcome

    monkeypatch.setattr(detect_cli, "MatpredictConfig", SimpleNamespace(from_env=lambda repo_root: SimpleNamespace(db_root=Path("db"))))
    seen: list = []

    def fake_run_pipeline(**kwargs):
        seen.append(kwargs.get("weak_flank_cap_policy", "MISSING"))
        return DetectionOutcome(results=[])

    monkeypatch.setattr(detect_cli, "build_reference_fasta", lambda db_root, out, family_keys=None, **kw: out)
    monkeypatch.setattr(detect_cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(detect_cli, "write_detection_gff3", lambda *a, **k: None)
    monkeypatch.setattr(detect_cli, "write_detection_report", lambda outcome, path, run=None: None)
    base = dict(genome="g.fa", proteins=None, taxid=None, out_dir=str(tmp_path / "out"), evidence_diagnostics=None, min_hits=1,
                min_identity=None, require_core_role=False, exclude_records="", emit_cds_fasta=False, phylum="Mucoromycota")
    assert detect_cli._cmd_detect(SimpleNamespace(**base)) == 0
    assert detect_cli._cmd_detect(SimpleNamespace(**base, weak_flank_cap_policy="on")) == 0
    assert seen == ["off", "on"]


def test_no_preset_turns_the_cap_on():
    # only the v0.9.1 preset turns this rule on (curator, 2026-11-07)
    assert {n for n, v in detect_cli.PRESETS.items() if "weak_flank_cap_policy" in v} == {"v0.9.1"}
    assert detect_cli.PRESETS["v0.9.1"]["weak_flank_cap_policy"] == "on"
