"""A genome whose phylum has no curated family is not searched by default.

Curator's ruling 2026-09-26. The chytrid negative control ran `exhaustive`
(every family of every phylum): 0 calls in 19 finished genomes, a median
1,866 s each, and 6 of 25 hit the 60-minute timeout. The default is now
`not_searched`: a report that says so and runs no search. `--exhaustive` (or
`--phylum`) still asks for a search explicitly.
"""
from pathlib import Path
from types import SimpleNamespace

import yaml

from MATPredict.detect.family_registry import RoutingDecision
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.report import write_detection_report
from MATPredict.detect.rollout_aggregate import aggregate_reports

from tests.detect.test_pipeline import _write_order, _write_record


def _boom(*a, **k):
    raise AssertionError("a not_searched run must not search or polish")


def test_run_pipeline_searches_nothing_when_not_searched(tmp_path):
    _write_order(tmp_path)
    _write_record(tmp_path)
    routing = RoutingDecision(families=[], routing_mode="not_searched",
                              phylum="Chytridiomycota")
    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=12345,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=_boom, search_localize=_boom,
        polish_with_exonerate=_boom, polish_with_miniprot=_boom,
        routing=routing, genetic_code_resolver=_boom,
    )
    assert outcome.results == []
    assert outcome.not_detected == []
    assert outcome.families_attempted == []
    assert outcome.routing_mode == "not_searched"
    assert "Chytridiomycota" in outcome.not_searched_reason


def test_the_report_states_not_searched(tmp_path):
    _write_order(tmp_path)
    _write_record(tmp_path)
    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=12345,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=_boom, search_localize=_boom,
        routing=RoutingDecision(families=[], routing_mode="not_searched"),
    )
    path = tmp_path / "report.yaml"
    write_detection_report(outcome, path)
    doc = yaml.safe_load(path.read_text())
    assert doc["routing_mode"] == "not_searched"
    assert doc["not_searched_reason"]
    assert doc["detected"] == [] and doc["families_attempted"] == []


def _cli(monkeypatch, tmp_path, routed_to, recorded):
    import MATPredict.detect.cli as detect_cli
    from MATPredict.detect.pipeline import DetectionOutcome

    monkeypatch.setattr(
        detect_cli, "MatpredictConfig",
        SimpleNamespace(from_env=lambda repo_root: SimpleNamespace(db_root=tmp_path / "db")),
    )
    monkeypatch.setattr(detect_cli, "load_all_families", lambda db_root: [])

    def fake_route(taxid, families, phylum=None, exhaustive=False):
        recorded["exhaustive"] = exhaustive
        return routed_to

    monkeypatch.setattr(detect_cli, "route", fake_route)

    def fake_build(*a, **k):
        recorded["built_reference"] = True
        return tmp_path / "ref.faa"

    monkeypatch.setattr(detect_cli, "build_reference_fasta", fake_build)
    monkeypatch.setattr(
        detect_cli, "run_pipeline",
        lambda **kw: DetectionOutcome(results=[], routing_mode=kw["routing"].routing_mode),
    )
    monkeypatch.setattr(detect_cli, "write_detection_gff3", lambda outcome, path, **k: None)
    return detect_cli


def _args(tmp_path, **extra):
    base = dict(
        genome="g.fa", proteins=None, taxid=12345, out_dir=str(tmp_path / "out"),
        evidence_diagnostics=None, min_hits=1, min_identity=None,
        require_core_role=False, exclude_records="", emit_cds_fasta=False, phylum=None,
    )
    base.update(extra)
    return SimpleNamespace(**base)


def test_the_cli_writes_a_not_searched_report_without_building_a_query_set(monkeypatch, tmp_path):
    recorded: dict = {}
    cli = _cli(monkeypatch, tmp_path,
               RoutingDecision(families=[], routing_mode="not_searched"), recorded)
    assert cli._cmd_detect(_args(tmp_path)) == 0
    assert "built_reference" not in recorded
    assert recorded["exhaustive"] is False
    doc = yaml.safe_load((tmp_path / "out" / "detection_report.yaml").read_text())
    assert doc["routing_mode"] == "not_searched"


def test_the_cli_passes_the_exhaustive_override(monkeypatch, tmp_path):
    recorded: dict = {}
    cli = _cli(monkeypatch, tmp_path,
               RoutingDecision(families=[], routing_mode="exhaustive"), recorded)
    cli._cmd_detect(_args(tmp_path, exhaustive=True))
    assert recorded["exhaustive"] is True


def test_the_exhaustive_flag_parses_and_defaults_off():
    import argparse
    from MATPredict.detect.cli import register_subcommands
    parser = argparse.ArgumentParser()
    register_subcommands(parser.add_subparsers(dest="command"))
    assert parser.parse_args(["detect"]).exhaustive is False
    assert parser.parse_args(["detect", "--exhaustive"]).exhaustive is True


def test_the_rollout_counts_not_searched_genomes_separately(tmp_path):
    def write(name, doc):
        d = tmp_path / name
        d.mkdir()
        p = d / "detection_report.yaml"
        p.write_text(yaml.safe_dump(doc))
        return p

    searched = write("111_GCA_1", {
        "routing_mode": "lineage", "families_attempted": ["P:aLocus"],
        "detected": [{"family": "P:aLocus", "confidence": "high"}], "not_detected": [],
    })
    skipped = write("111_GCA_2", {
        "routing_mode": "not_searched", "not_searched_reason": "no curated family",
        "families_attempted": [], "detected": [], "not_detected": [],
    })
    summary = aggregate_reports([searched, skipped], lineage_resolver=lambda _t: "SameOrder")
    assert summary.not_searched == ["111_GCA_2"]
    assert summary.genome_errors == []
    # A genome that was never searched cannot be an anomaly for missing a family.
    assert summary.anomalies == []
    assert summary.to_doc()["not_searched"] == ["111_GCA_2"]
