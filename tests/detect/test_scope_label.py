"""A Basidiomycota MAT call outside its scope on weak genes only is labelled unverified (`scope_label`, #170)."""
import argparse
import inspect
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from MATPredict.detect import cli as detect_cli
from MATPredict.detect import scope_label as sl
from MATPredict.detect.family_registry import Family, FamilyKey, _gene_alias_map
from MATPredict.detect.pipeline import DetectionResult, run_pipeline

MAT, HD = FamilyKey("Basidiomycota", "MAT"), FamilyKey("Basidiomycota", "HD")
TREMELLO, OTHER = 5580, 5338


def _family(key, scope):
    genes = [{"name": "X", "role": "core_MAT"}]
    return Family(key=key, vocabulary_type="enum", idiomorph_values=["a", "alpha"], idiomorph_pattern=None, genes=genes,
                  taxonomic_scope=scope, gene_aliases=_gene_alias_map(genes))


FAMILIES = [_family(MAT, [TREMELLO]), _family(HD, [OTHER])]


def _call(key, genes, confidence="high"):
    return DetectionResult(family_key=key, contig="c1", start=1, end=100, confidence=confidence, idiomorph="alpha", ambiguous_with=[],
                           genes_found=list(genes), genes_missing=[], fragmented=False, gene_evidence=[])


def test_defaults_are_off():
    assert sl.DEFAULT_OUT_OF_SCOPE_LABEL_POLICY == "off"
    assert inspect.signature(run_pipeline).parameters["out_of_scope_label_policy"].default == "off"
    assert detect_cli._out_of_scope_label_policy_from_args(argparse.Namespace()) == "off"


def test_unknown_policy_is_an_error():
    with pytest.raises(ValueError):
        sl.check_policy("yes")


def test_off_returns_the_calls_untouched():
    r = _call(MAT, ["STE3", "PAN6"])
    assert sl.apply([r], FAMILIES, {OTHER}, "off") == [r]


@pytest.mark.parametrize("genes", [["STE3", "PAN6"], ["STE3", "RPL39"], ["STE3"], ["PAN6", "RPL39", "STE3"]])
def test_weak_genes_only_outside_the_scope_is_unconfirmed_and_nothing_else_changes(genes):
    r = _call(MAT, genes)
    (out,) = sl.apply([r], FAMILIES, {OTHER}, "on")
    assert out.family_status["status"] == "family_unconfirmed" and out.family_status["reason"] == sl.REASON
    assert (out.idiomorph, out.confidence, out.start, out.end, out.genes_found) == (r.idiomorph, r.confidence, r.start, r.end, r.genes_found)


@pytest.mark.parametrize("genes,lineage,key", [
    (["STE3", "PAN6"], {TREMELLO}, MAT),                 # inside the scope
    (["STE3", "PAN6"], set(), MAT),                      # unknown lineage: never flagged
    (["STE3", "PAN6", "SXI1"], {OTHER}, MAT),            # a gene beyond the weak ones
    (["MFalpha", "PAN6", "STE3"], {OTHER}, MAT),
    (["FAO1", "RPL39"], {OTHER}, MAT),
    ([], {OTHER}, MAT),
    (["STE3", "PAN6"], {TREMELLO}, HD),                  # another family is not touched
])
def test_everything_else_is_untouched(genes, lineage, key):
    r = _call(key, genes)
    assert sl.apply([r], FAMILIES, lineage, "on") == [r]


def test_an_existing_family_status_is_kept():
    r = replace(_call(MAT, ["STE3"]), family_status={"status": "family_confirmed"})
    assert sl.apply([r], FAMILIES, {OTHER}, "on") == [r]


def test_cli_default_and_passthrough(monkeypatch, tmp_path):
    from MATPredict.detect.pipeline import DetectionOutcome

    monkeypatch.setattr(detect_cli, "MatpredictConfig", SimpleNamespace(from_env=lambda repo_root: SimpleNamespace(db_root=Path("db"))))
    seen: list = []

    def fake_run_pipeline(**kwargs):
        seen.append(kwargs.get("out_of_scope_label_policy", "MISSING"))
        return DetectionOutcome(results=[])

    monkeypatch.setattr(detect_cli, "build_reference_fasta", lambda db_root, out, family_keys=None, **kw: out)
    monkeypatch.setattr(detect_cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(detect_cli, "write_detection_gff3", lambda *a, **k: None)
    monkeypatch.setattr(detect_cli, "write_detection_report", lambda outcome, path, run=None: None)
    base = dict(genome="g.fa", proteins=None, taxid=None, out_dir=str(tmp_path / "out"), evidence_diagnostics=None, min_hits=1,
                min_identity=None, require_core_role=False, exclude_records="", emit_cds_fasta=False, phylum="Mucoromycota")
    assert detect_cli._cmd_detect(SimpleNamespace(**base)) == 0
    assert detect_cli._cmd_detect(SimpleNamespace(**base, out_of_scope_label_policy="on")) == 0
    assert seen == ["off", "on"]


def test_no_preset_sets_the_policy():
    assert all("out_of_scope_label_policy" not in v for v in detect_cli.PRESETS.values())
