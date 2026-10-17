"""The confidence cap (`confidence_cap`, `--confidence-cap-policy`).

Pure-function tests of `confidence_cap.evaluate`, then the same scenarios through
`run_pipeline` with stubbed polishing: a two-gene, flank-less locus (a MATsc-like
`idiomorph_gene_only` call, `high` when both genes are modelled).
"""
from __future__ import annotations

import argparse
import inspect

import pytest

from MATPredict.detect import confidence_cap as cc
from MATPredict.detect import pipeline as pipeline_module
from MATPredict.detect.cli import _confidence_cap_policy_from_args
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.polish import ExonSpan, PolishModel
from MATPredict.detect.search import SearchHit


def S(gene="G1", cov=100.0, ident=90.0, bits=300.0, edge=50000):
    return cc.CapStat(gene, cov, ident, bits, edge)


def test_default_off_and_check_policy():
    assert cc.DEFAULT_CONFIDENCE_CAP_POLICY == "off"
    assert _confidence_cap_policy_from_args(argparse.Namespace()) == "off"
    assert inspect.signature(run_pipeline).parameters["confidence_cap_policy"].default == "off"
    for p in cc.CONFIDENCE_CAP_POLICIES:
        assert cc.check_policy(p) == p
    with pytest.raises(ValueError):
        cc.check_policy("bogus")


def test_off_never_caps():
    assert cc.evaluate("high", [S(ident=30.0)], 0, policy="off") is None


def test_only_high_is_capped():
    for tier in ("medium", "low"):
        assert cc.evaluate(tier, [S(ident=30.0)], 0) is None


def test_weak_identity_without_flank_caps():
    d = cc.evaluate("high", [S(ident=39.06)], 0)
    assert d.criteria == (cc.CRITERION_WEAK,)
    assert d.note["to"] == "medium" and d.note["best_core_identity"] == 39.06 and d.note["flank_genes"] == 0


def test_weak_identity_with_a_flank_gene_does_not_cap():
    assert cc.evaluate("high", [S(ident=39.06)], 1) is None


def test_identity_at_or_above_threshold_does_not_cap():
    assert cc.evaluate("high", [S(ident=45.0)], 0) is None
    assert cc.evaluate("high", [S(ident=99.0)], 0) is None


def test_best_core_is_the_highest_bitscore_not_the_highest_identity():
    stats = [S("A", ident=40.0, bits=500.0), S("B", ident=95.0, bits=100.0)]
    assert cc.evaluate("high", stats, 0).note["best_core_gene"] == "A"
    stats = [S("A", ident=40.0, bits=100.0), S("B", ident=95.0, bits=500.0)]
    assert cc.evaluate("high", stats, 0) is None


def test_unknown_identity_or_no_core_never_caps():
    assert cc.evaluate("high", [S(ident=None)], 0) is None
    assert cc.evaluate("high", [], 0) is None


def test_partial_at_the_contig_edge_caps():
    d = cc.evaluate("high", [S(cov=60.0, edge=120)], 1)
    assert d.criteria == (cc.CRITERION_PARTIAL,)


def test_partial_away_from_the_contig_edge_does_not_cap():
    # a domain-only model against a long reference is ordinary (49% of v0.6.0 high calls)
    assert cc.evaluate("high", [S(cov=30.0, edge=50000)], 1) is None
    assert cc.evaluate("high", [S(cov=30.0, edge=None)], 1) is None


def test_edge_partial_with_a_full_length_core_gene_does_not_cap():
    stats = [S("A", cov=60.0, edge=120, bits=400.0), S("B", cov=95.0, edge=90000, bits=100.0)]
    assert cc.evaluate("high", stats, 1) is None


def test_full_or_80_percent_is_not_partial():
    assert cc.evaluate("high", [S(cov=80.0, edge=10)], 1) is None
    assert cc.evaluate("high", [S(cov=79.9, edge=10)], 1) is not None


def test_both_criteria_are_listed():
    d = cc.evaluate("high", [S(cov=50.0, ident=40.0, edge=10)], 0)
    assert d.criteria == (cc.CRITERION_PARTIAL, cc.CRITERION_WEAK)


def test_weak_only_policy_ignores_the_partial_criterion():
    assert cc.evaluate("high", [S(cov=50.0, edge=10)], 1, policy="weak_only") is None
    d = cc.evaluate("high", [S(cov=50.0, ident=40.0, edge=10)], 0, policy="weak_only")
    assert d.criteria == (cc.CRITERION_WEAK,)


# ---- through run_pipeline ---------------------------------------------------

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


def _run(tmp_path, monkeypatch, *, policy, ident=39.0, a1_start=30000, contig_len=60000, aa=290):
    _write_records(tmp_path)
    (tmp_path / "genome.fa").write_text(">c1\n" + "A" * contig_len + "\n")
    monkeypatch.setattr(pipeline_module, "_translate_model", lambda *a, **k: "MKV")
    starts = {"ALPHA1": (a1_start, "r1"), "ALPHA2": (a1_start + 3000, "r2")}

    def model(gene_name, method):
        start, rid = starts[gene_name]
        return PolishModel(gene_name, KEY, "core_MAT", "c1", start, start + aa * 3 - 1, "+",
                           [ExonSpan(start, start + aa * 3 - 1)], ident, rid, method)

    hits = [SearchHit(KEY, g, "core_MAT", "c1", s, s + aa * 3 - 1, "+", ident, rid, "tblastn_genome", bitscore=200.0)
            for g, (s, rid) in starts.items()]
    out = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None, db_root=tmp_path,
        reference_fasta=tmp_path / "reference.faa", search_localize=lambda *a, **k: hits,
        polish_with_exonerate=lambda *, gene_name, window, **kw: model(gene_name, "exonerate_refine"),
        polish_with_miniprot=lambda *, gene_name, window, **kw: model(gene_name, "miniprot_refine"),
        confidence_cap_policy=policy,
    )
    assert len(out.results) == 1
    return out.results[0]


def test_baseline_two_gene_flankless_locus_is_high(tmp_path, monkeypatch):
    r = _run(tmp_path, monkeypatch, policy="off")
    assert (r.locus_class, r.confidence, r.idiomorph) == ("idiomorph_gene_only", "high", "alpha")
    assert r.confidence_cap is None


def test_weak_gene_only_call_drops_to_medium_and_nothing_else_changes(tmp_path, monkeypatch):
    off = _run(tmp_path / "a" if (tmp_path / "a").mkdir() is None else None, monkeypatch, policy="off", ident=39.0)
    on = _run(tmp_path / "b" if (tmp_path / "b").mkdir() is None else None, monkeypatch, policy="partial_weak", ident=39.0)
    assert (off.confidence, on.confidence) == ("high", "medium")
    assert on.confidence_cap["criteria"] == ["weak_identity"] and on.confidence_cap["from"] == "high"
    for f in ("idiomorph", "locus_class", "genes_found", "genes_missing", "start", "end", "contig", "polished_genes",
              "idiomorph_candidates", "idiomorph_margin", "detection_pass"):
        assert getattr(on, f) == getattr(off, f), f


def test_strong_gene_only_call_stays_high_under_on(tmp_path, monkeypatch):
    r = _run(tmp_path, monkeypatch, policy="partial_weak", ident=98.0)
    assert r.confidence == "high" and r.confidence_cap is None


def test_report_carries_the_note_only_when_set(tmp_path, monkeypatch):
    from MATPredict.detect.report import _result_doc
    (tmp_path / "a").mkdir(); (tmp_path / "b").mkdir()
    off = _result_doc(_run(tmp_path / "a", monkeypatch, policy="off", ident=39.0))
    on = _result_doc(_run(tmp_path / "b", monkeypatch, policy="partial_weak", ident=39.0))
    assert "confidence_cap" not in off and on["confidence_cap"]["to"] == "medium"
