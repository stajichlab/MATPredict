"""The weak-fragment tier cap (`weak_fragment_cap`, `--weak-fragment-cap-policy`, #275).

Pure-function tests of `weak_fragment_cap.evaluate`, then the same scenarios through `run_pipeline` with stubbed polishing (the scope
constants are pointed at the test family), then the CLI default and passthrough.
"""
from __future__ import annotations

import argparse
import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

from MATPredict.detect import pipeline as pipeline_module
from MATPredict.detect import weak_fragment_cap as wfc
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.polish import ExonSpan, PolishModel
from MATPredict.detect.search import SearchHit

FAM = ("Ascomycota", "MAT")


def C(gene="MAT1-2-1", status="polished_disagree", bits=45.0):
    return wfc.CoreStat(gene, status, bits)


def test_default_is_off_and_check_policy():
    assert wfc.DEFAULT_WEAK_FRAGMENT_CAP_POLICY == "off"
    assert inspect.signature(run_pipeline).parameters["weak_fragment_cap_policy"].default == "off"
    for p in wfc.WEAK_FRAGMENT_CAP_POLICIES:
        assert wfc.check_policy(p) == p
    with pytest.raises(ValueError):
        wfc.check_policy("bogus")


def test_off_never_caps():
    assert wfc.evaluate("high", "MAT1-2", FAM, [C()], policy="off") is None


@pytest.mark.parametrize("status", ["polished_disagree", "polished_single", "unpolished"])
@pytest.mark.parametrize("tier", ["medium", "high"])
def test_a_weak_strongest_core_gene_caps_medium_and_high_to_low(status, tier):
    d = wfc.evaluate(tier, "MAT1-2", FAM, [C(status=status, bits=45.0)])
    assert d.note["from"] == tier and d.note["to"] == "low"
    assert d.note["strongest_core_gene"] == "MAT1-2-1" and d.note["strongest_core_bitscore"] == 45.0 and d.note["threshold_bitscore"] == 50.0


def test_low_is_untouched():
    assert wfc.evaluate("low", "MAT1-2", FAM, [C()]) is None


def test_a_both_tools_agree_model_is_never_a_weak_fragment():
    assert wfc.evaluate("high", "MAT1-2", FAM, [C(status="polished_agree", bits=30.0)]) is None


def test_fifty_bits_or_more_is_not_weak():
    assert wfc.evaluate("high", "MAT1-2", FAM, [C(bits=50.0)]) is None
    assert wfc.evaluate("high", "MAT1-2", FAM, [C(bits=49.9)]) is not None


def test_the_strongest_core_gene_decides_not_the_weakest():
    stats = [C("MAT1-1-3", "polished_single", 30.0), C("MAT1-2-1", "polished_agree", 210.0)]
    assert wfc.evaluate("high", "MAT1-2", FAM, stats) is None            # an agreeing, strong gene is the strongest
    stats = [C("MAT1-1-3", "unpolished", 90.0), C("MAT1-2-1", "polished_single", 44.0)]
    d = wfc.evaluate("high", "MAT1-2", FAM, stats)                        # polished beats unpolished first, then the bitscore
    assert d.note["strongest_core_gene"] == "MAT1-2-1"


def test_flank_genes_are_not_core_genes():
    stats = [wfc.CoreStat("SLA2", "polished_agree", 1700.0), C("MAT1-2-1", "unpolished", 40.0)]
    assert wfc.evaluate("high", "MAT1-2", FAM, stats) is not None


def test_only_the_ascomycota_mat_family_and_called_idiomorphs_are_in_scope():
    assert wfc.evaluate("high", "MAT1-2", ("Ascomycota", "MATsc"), [C()]) is None
    assert wfc.evaluate("high", "MAT1-2", ("Basidiomycota", "MAT"), [C()]) is None
    assert wfc.evaluate("high", "undetermined", FAM, [C()]) is None
    assert wfc.evaluate("high", None, FAM, [C()]) is None
    assert wfc.evaluate("high", "MAT1-2", FAM, []) is None


# ---- through run_pipeline ---------------------------------------------------------------------------------------------------------

ORDER = (
    "phylum: P\nloci:\n  - locus_name: MATX\n    vocabulary_type: enum\n"
    "    idiomorph_values: [alpha, a]\n    taxonomic_scope: [1]\n    genes:\n"
    "      - {name: ALPHA1, role: core_MAT, present_in_idiomorphs: [\"alpha\"]}\n"
    "      - {name: ALPHA2, role: core_MAT, present_in_idiomorphs: [\"alpha\"]}\n"
    "      - {name: A1, role: core_MAT, present_in_idiomorphs: [\"a\"]}\n"
)
KEY = pipeline_module.FamilyKey("P", "MATX")


def _write_records(tmp_path):
    (tmp_path / "P").mkdir(exist_ok=True)
    (tmp_path / "P" / "order.yml").write_text(ORDER)
    for rid, gene, idiom in (("r1", "ALPHA1", "alpha"), ("r2", "ALPHA2", "alpha"), ("r3", "A1", "a")):
        d = tmp_path / "P" / "Fam" / rid
        d.mkdir(parents=True, exist_ok=True)
        (d / "metadata.yaml").write_text(f"record_id: {rid}\nmating_type: {{locus_name: MATX, idiomorphs: [{idiom}]}}\n")
        (d / "proteins.faa").write_text(f">{rid}|gene_index=0|name={gene}|role=core_MAT\n" + "M" * 300 + "\n")


def _run(tmp_path, monkeypatch, *, policy, bits=200.0, both_tools=True, aa=290, ident=90.0):
    monkeypatch.setattr(wfc, "SCOPE_FAMILY", ("P", "MATX"))
    monkeypatch.setattr(wfc, "CALLED_IDIOMORPHS", ("alpha", "a"))
    monkeypatch.setattr(wfc, "CORE_GENE_PREFIX", "ALPHA")
    _write_records(tmp_path)
    (tmp_path / "genome.fa").write_text(">c1\n" + "A" * 60000 + "\n")
    monkeypatch.setattr(pipeline_module, "_translate_model", lambda *a, **k: "MKV")
    starts = {"ALPHA1": (30000, "r1"), "ALPHA2": (33000, "r2")}

    def model(gene_name, method):
        start, rid = starts[gene_name]
        return PolishModel(gene_name, KEY, "core_MAT", "c1", start, start + aa * 3 - 1, "+",
                           [ExonSpan(start, start + aa * 3 - 1)], ident, rid, method)

    hits = [SearchHit(KEY, g, "core_MAT", "c1", s, s + aa * 3 - 1, "+", ident, rid, "tblastn_genome", bitscore=bits)
            for g, (s, rid) in starts.items()]
    out = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None, db_root=tmp_path,
        reference_fasta=tmp_path / "reference.faa", search_localize=lambda *a, **k: hits,
        polish_with_exonerate=lambda *, gene_name, window, **kw: model(gene_name, "exonerate_refine"),
        polish_with_miniprot=(lambda *, gene_name, window, **kw: model(gene_name, "miniprot_refine")) if both_tools
        else (lambda **kw: None),
        weak_fragment_cap_policy=policy,
    )
    assert len(out.results) == 1
    return out.results[0]


def test_a_single_tool_weak_call_drops_to_low_and_nothing_else_changes(tmp_path, monkeypatch):
    (tmp_path / "a").mkdir(); (tmp_path / "b").mkdir()
    off = _run(tmp_path / "a", monkeypatch, policy="off", bits=40.0, both_tools=False)
    on = _run(tmp_path / "b", monkeypatch, policy="on", bits=40.0, both_tools=False)
    assert off.confidence in ("medium", "high") and on.confidence == "low" and on.weak_fragment_cap["from"] == off.confidence
    assert off.weak_fragment_cap is None
    for f in ("idiomorph", "locus_class", "genes_found", "genes_missing", "start", "end", "contig", "polished_genes",
              "idiomorph_candidates", "idiomorph_margin", "detection_pass"):
        assert getattr(on, f) == getattr(off, f), f


def test_a_strong_call_and_a_both_tools_agree_call_keep_their_tier(tmp_path, monkeypatch):
    (tmp_path / "a").mkdir(); (tmp_path / "b").mkdir()
    strong = _run(tmp_path / "a", monkeypatch, policy="on", bits=200.0, both_tools=False)
    agree = _run(tmp_path / "b", monkeypatch, policy="on", bits=40.0, both_tools=True)
    assert strong.weak_fragment_cap is None and agree.weak_fragment_cap is None
    assert strong.confidence != "low" and agree.confidence != "low"


def test_report_carries_the_note_only_when_set(tmp_path, monkeypatch):
    from MATPredict.detect.report import _result_doc
    (tmp_path / "a").mkdir(); (tmp_path / "b").mkdir()
    off = _result_doc(_run(tmp_path / "a", monkeypatch, policy="off", bits=40.0, both_tools=False))
    on = _result_doc(_run(tmp_path / "b", monkeypatch, policy="on", bits=40.0, both_tools=False))
    assert "weak_fragment_cap" not in off and on["weak_fragment_cap"]["to"] == "low"


# ---- CLI ---------------------------------------------------------------------------------------------------------------------------

def test_cli_default_and_passthrough(monkeypatch, tmp_path):
    from MATPredict.detect import cli as detect_cli
    from MATPredict.detect.pipeline import DetectionOutcome

    assert detect_cli._weak_fragment_cap_policy_from_args(argparse.Namespace()) == "off"
    monkeypatch.setattr(detect_cli, "MatpredictConfig", SimpleNamespace(from_env=lambda repo_root: SimpleNamespace(db_root=Path("db"))))
    seen: list = []

    def fake_run_pipeline(**kwargs):
        seen.append(kwargs.get("weak_fragment_cap_policy", "MISSING"))
        return DetectionOutcome(results=[])

    monkeypatch.setattr(detect_cli, "build_reference_fasta", lambda db_root, out, family_keys=None, **kw: out)
    monkeypatch.setattr(detect_cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(detect_cli, "write_detection_gff3", lambda *a, **k: None)
    monkeypatch.setattr(detect_cli, "write_detection_report", lambda outcome, path, run=None: None)
    base = dict(genome="g.fa", proteins=None, taxid=None, out_dir=str(tmp_path / "out"), evidence_diagnostics=None, min_hits=1,
                min_identity=None, require_core_role=False, exclude_records="", emit_cds_fasta=False, phylum="Mucoromycota")
    assert detect_cli._cmd_detect(SimpleNamespace(**base)) == 0
    assert detect_cli._cmd_detect(SimpleNamespace(**base, weak_fragment_cap_policy="on")) == 0
    assert seen == ["off", "on"]


def test_no_preset_turns_the_cap_on():
    from MATPredict.detect.cli import PRESETS
    assert all("weak_fragment_cap_policy" not in v for v in PRESETS.values())
