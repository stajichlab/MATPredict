"""`--withheld-detail-policy`: a report-only per-gene summary on withheld loci (the one-polished-gene measurement; `why no call` detail)."""
import argparse
import inspect

import yaml

from MATPredict.detect import cli as detect_cli
from MATPredict.detect import report as R
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import DetectionOutcome, DetectionResult, GeneEvidence

KEY = FamilyKey("Ascomycota", "MAT")


def _gene(name, role, status, identity=55.123, coverage=80.0, bitscore=120.04, lost_to=None):
    return GeneEvidence(gene_name=name, role=role, contig="c1", start=100, end=400, strand="+", identity=identity, coverage=coverage,
                        reference_record_id="r", method="tblastn_genome", status=status, bitscore=bitscore, lost_to=lost_to)


def _suppressed():
    genes = [_gene("MAT1-1-1", "core_MAT", "polished_agree"), _gene("SLA2", "flanking_conserved", "unpolished", identity=None),
             _gene("MAT1-2-1", "core_MAT", "polished_agree", lost_to="MAT1-1-1")]
    return DetectionResult(family_key=KEY, contig="c1", start=1, end=900, confidence="low", idiomorph="MAT1-1", ambiguous_with=[],
                           genes_found=["MAT1-1-1", "SLA2"], genes_missing=[], fragmented=False, gene_evidence=genes, polished_genes=1,
                           withheld_reason="modelled_gene_bar")


def _report(tmp_path, policy=None):
    out = tmp_path / f"r_{policy}.yaml"
    run = None if policy is None else {"parameters": {"withheld_detail_policy": policy}}
    R.write_detection_report(DetectionOutcome(results=[], suppressed_loci=[_suppressed()]), out, run=run)
    return yaml.safe_load(out.read_text())["suppressed_loci"][0]


def test_defaults_are_off():
    assert R.DEFAULT_WITHHELD_DETAIL_POLICY == "off"
    assert detect_cli._withheld_detail_policy_from_args(argparse.Namespace()) == "off"


def test_no_genes_key_by_default_or_when_off(tmp_path):
    assert "genes" not in _report(tmp_path) and "genes" not in _report(tmp_path, "off")


def test_on_lists_the_live_genes_with_role_status_and_numbers(tmp_path):
    entry = _report(tmp_path, "on")
    assert entry["genes"] == [
        {"gene": "MAT1-1-1", "role": "core_MAT", "status": "polished_agree", "identity": 55.12, "coverage": 80.0, "bitscore": 120.0},
        {"gene": "SLA2", "role": "flanking_conserved", "status": "unpolished", "identity": None, "coverage": 80.0, "bitscore": 120.0},
    ]


def test_on_changes_nothing_else_in_the_entry(tmp_path):
    off, on = _report(tmp_path, "off"), _report(tmp_path, "on")
    assert {k: v for k, v in on.items() if k != "genes"} == off


def test_the_policy_is_a_cli_flag_default_off_and_in_no_preset():
    parser = argparse.ArgumentParser()
    detect_cli.register_subcommands(parser.add_subparsers(dest="cmd"))
    args = parser.parse_args(["detect", "--genome", "g.fa", "--taxid", "1"])
    assert args.withheld_detail_policy == "off"
    assert parser.parse_args(["detect", "--genome", "g.fa", "--taxid", "1", "--withheld-detail-policy", "on"]).withheld_detail_policy == "on"
    assert all("withheld_detail_policy" not in v for v in detect_cli.PRESETS.values())
