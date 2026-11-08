"""The single-region tier cap (`single_region_cap`, `--single-region-cap-policy`, #167)."""
import argparse
import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

from MATPredict.detect import cli as detect_cli
from MATPredict.detect import single_region_cap as sr
from MATPredict.detect.pipeline import run_pipeline


def G(gene, role="core_MAT", status="polished_single", contig="c1", start=1000, end=1900):
    return sr.GeneStat(gene, role, status, contig, start, end)


def test_defaults_are_off():
    assert sr.DEFAULT_SINGLE_REGION_CAP_POLICY == "off"
    assert inspect.signature(run_pipeline).parameters["single_region_cap_policy"].default == "off"
    assert detect_cli._single_region_cap_policy_from_args(argparse.Namespace()) == "off"


def test_unknown_policy_is_an_error():
    with pytest.raises(ValueError):
        sr.check_policy("always")


def test_overlapping_hd1_and_hd2_are_one_region():
    assert sr.count_regions([G("HD1", start=1000, end=1900), G("HD2", start=1100, end=1950)]) == 1


def test_separate_models_are_two_regions():
    assert sr.count_regions([G("HD1", start=1000, end=1900), G("HD2", start=5000, end=5900)]) == 2


def test_different_contigs_are_two_regions_and_a_small_overlap_is_not_one():
    assert sr.count_regions([G("HD1", contig="a"), G("HD2", contig="b")]) == 2
    assert sr.count_regions([G("HD1", start=1000, end=1900), G("HD2", start=1800, end=2700)]) == 2   # 101 of 901 bp


def test_off_never_caps():
    assert sr.evaluate("high", "Basidiomycota", [G("HD1"), G("HD2", start=1100, end=1950)], policy="off") is None


def test_one_shared_region_without_a_polished_flank_caps_a_high_call():
    note = sr.evaluate("high", "Basidiomycota", [G("HD1"), G("HD2", start=1100, end=1950)], policy="on")
    assert note["from"] == "high" and note["to"] == "medium" and note["distinct_core_regions"] == 1 and note["core_genes"] == ["HD1", "HD2"]


def test_a_single_gene_call_is_capped():
    assert sr.evaluate("high", "Basidiomycota", [G("bE", end=1100)], policy="on")["distinct_core_regions"] == 1


@pytest.mark.parametrize("stats,confidence,phylum", [
    ([G("HD1"), G("HD2", start=5000, end=5900)], "high", "Basidiomycota"),                                   # two distinct regions
    ([G("HD1"), G("PAP", role="flanking_conserved", status="polished_agree", start=9000, end=9500)], "high", "Basidiomycota"),  # polished flank
    ([G("HD1")], "medium", "Basidiomycota"),                                                                    # only high is capped
    ([G("HD1")], "high", "Ascomycota"),                                                                         # other phyla untouched
    ([G("PAP", role="flanking_conserved")], "high", "Basidiomycota"),                                           # no core gene
])
def test_everything_else_keeps_its_tier(stats, confidence, phylum):
    assert sr.evaluate(confidence, phylum, stats, policy="on") is None


def test_an_unpolished_flank_gene_does_not_save_the_call():
    flank = G("PAP", role="flanking_conserved", status="unpolished", start=9000, end=9500)
    assert sr.evaluate("high", "Basidiomycota", [G("HD1"), flank], policy="on") is not None


def test_cli_default_and_passthrough(monkeypatch, tmp_path):
    from MATPredict.detect.pipeline import DetectionOutcome

    monkeypatch.setattr(detect_cli, "MatpredictConfig", SimpleNamespace(from_env=lambda repo_root: SimpleNamespace(db_root=Path("db"))))
    seen: list = []

    def fake_run_pipeline(**kwargs):
        seen.append(kwargs.get("single_region_cap_policy", "MISSING"))
        return DetectionOutcome(results=[])

    monkeypatch.setattr(detect_cli, "build_reference_fasta", lambda db_root, out, family_keys=None, **kw: out)
    monkeypatch.setattr(detect_cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(detect_cli, "write_detection_gff3", lambda *a, **k: None)
    monkeypatch.setattr(detect_cli, "write_detection_report", lambda outcome, path, run=None: None)
    base = dict(genome="g.fa", proteins=None, taxid=None, out_dir=str(tmp_path / "out"), evidence_diagnostics=None, min_hits=1,
                min_identity=None, require_core_role=False, exclude_records="", emit_cds_fasta=False, phylum="Mucoromycota")
    assert detect_cli._cmd_detect(SimpleNamespace(**base)) == 0
    assert detect_cli._cmd_detect(SimpleNamespace(**base, single_region_cap_policy="on")) == 0
    assert seen == ["off", "on"]


def test_no_preset_turns_the_cap_on():
    assert all("single_region_cap_policy" not in v for v in detect_cli.PRESETS.values())
