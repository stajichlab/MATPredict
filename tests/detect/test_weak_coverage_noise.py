"""An unpolished gene covering little of its reference is noise (`weak_coverage_noise`, `--weak-coverage-noise-policy`, #149)."""
import argparse
import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

from MATPredict.detect import cli as detect_cli
from MATPredict.detect import weak_coverage_noise as wc
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.search import SearchHit

K = FamilyKey("Ascomycota", "MAT")
OTHER = FamilyKey("Ascomycota", "MATsc")


def H(gene, coverage, key=K, superseded_by=None, bitscore=44.0):
    return SearchHit(family_key=key, gene_name=gene, role="core_MAT", contig="c1", start=1, end=100, strand="+", identity=40.0,
                     reference_record_id="r", method="tblastn_genome", coverage=coverage, bitscore=bitscore, superseded_by=superseded_by)


def test_defaults_are_off():
    assert wc.DEFAULT_WEAK_COVERAGE_NOISE_POLICY == "off"
    assert inspect.signature(run_pipeline).parameters["weak_coverage_noise_policy"].default == "off"
    assert detect_cli._weak_coverage_noise_policy_from_args(argparse.Namespace()) == "off"


def test_unknown_policy_is_an_error():
    with pytest.raises(ValueError):
        wc.check_policy("sometimes")


def test_a_gene_whose_every_hit_is_low_coverage_is_noise_even_above_the_bit_floor():
    assert wc.low_coverage_genes([H("APN2", 12.0, bitscore=44.0), H("APN2", 25.0, bitscore=42.0)], K) == {"APN2"}


def test_one_good_hit_keeps_the_gene():
    assert wc.low_coverage_genes([H("APN2", 12.0), H("APN2", 60.0)], K) == set()


def test_exactly_at_the_line_is_kept():
    assert wc.low_coverage_genes([H("APN2", 30.0)], K) == set()
    assert wc.low_coverage_genes([H("APN2", 29.9)], K) == {"APN2"}


def test_unknown_coverage_never_qualifies():
    assert wc.low_coverage_genes([H("APN2", None)], K) == set()
    assert wc.low_coverage_genes([H("APN2", 10.0), H("APN2", None)], K) == set()


def test_other_families_and_superseded_hits_are_ignored():
    assert wc.low_coverage_genes([H("APN2", 10.0, key=OTHER)], K) == set()
    assert wc.low_coverage_genes([H("APN2", 10.0, superseded_by="MAT1-1-1")], K) == set()


def test_cli_default_and_passthrough(monkeypatch, tmp_path):
    from MATPredict.detect.pipeline import DetectionOutcome

    monkeypatch.setattr(detect_cli, "MatpredictConfig", SimpleNamespace(from_env=lambda repo_root: SimpleNamespace(db_root=Path("db"))))
    seen: list = []

    def fake_run_pipeline(**kwargs):
        seen.append(kwargs.get("weak_coverage_noise_policy", "MISSING"))
        return DetectionOutcome(results=[])

    monkeypatch.setattr(detect_cli, "build_reference_fasta", lambda db_root, out, family_keys=None, **kw: out)
    monkeypatch.setattr(detect_cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(detect_cli, "write_detection_gff3", lambda *a, **k: None)
    monkeypatch.setattr(detect_cli, "write_detection_report", lambda outcome, path, run=None: None)
    base = dict(genome="g.fa", proteins=None, taxid=None, out_dir=str(tmp_path / "out"), evidence_diagnostics=None, min_hits=1,
                min_identity=None, require_core_role=False, exclude_records="", emit_cds_fasta=False, phylum="Mucoromycota")
    assert detect_cli._cmd_detect(SimpleNamespace(**base)) == 0
    assert detect_cli._cmd_detect(SimpleNamespace(**base, weak_coverage_noise_policy="on")) == 0
    assert seen == ["off", "on"]


def test_no_preset_sets_the_policy():
    assert all("weak_coverage_noise_policy" not in v for v in detect_cli.PRESETS.values())
