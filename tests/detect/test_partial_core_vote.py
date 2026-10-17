"""The partial-coverage vote rule (`partial_core`, `--partial-core-policy`).

Botrytis pattern (17 BioProject PRJNA624742 genomes stated MAT1-2, called
MAT1-1 high): a 214-aa 5'-truncated MAT1-1-1 remnant (61% of its 353-aa
reference, bitscore 416) outscores the real full-length 3-exon MAT1-2-1 (376
aa, 100% of its reference, bitscore 326) in the best-bitscore vote.

Pure-function tests of `partial_core.evaluate`, then the same scenarios through
`run_pipeline` (profile scoring stubbed: pyhmmer is not needed).
"""
from __future__ import annotations

import argparse

import pytest

from MATPredict.detect import partial_core as pc
from MATPredict.detect import pipeline as pipeline_module
from MATPredict.detect.cli import _partial_core_policy_from_args
from MATPredict.detect.pipeline import model_reference_coverage, run_pipeline
from MATPredict.detect.polish import ExonSpan, PolishModel
from MATPredict.detect.search import SearchHit

VOTE = [{"idiomorph": "MAT1-1", "score": 416.0}, {"idiomorph": "MAT1-2", "score": 326.0}]


def _stat(gene, idiom, cov, bits, edge=5000):
    return pc.CoreStat(gene, idiom, cov, bits, edge)


BOTRYTIS = [_stat("MAT1-1-1", "MAT1-1", 60.6, 416.0), _stat("MAT1-2-1", "MAT1-2", 100.0, 326.0)]


def _profiles(a, h):
    return lambda genes: {"MAT1-1-1": a, "MAT1-2-1": h}


def test_policy_constants_and_default_off():
    assert pc.DEFAULT_PARTIAL_CORE_POLICY == "off"
    assert _partial_core_policy_from_args(argparse.Namespace()) == "off"
    with pytest.raises(ValueError):
        pc.check_policy("bogus")


def test_botrytis_flips_when_profiles_agree():
    d = pc.evaluate(VOTE, BOTRYTIS, _profiles(125.0, 221.0))
    assert d.decision == pc.DECISION_FLIPPED and d.idiomorph == "MAT1-2"
    assert d.note["profile"]["margin_bits"] == 96.0
    assert d.note["winner_gene"] == "MAT1-1-1" and d.note["other_gene"] == "MAT1-2-1"


def test_botrytis_undetermined_when_profiles_disagree():
    # the profiles favour the partial MAT1-1-1 by more than the margin
    d = pc.evaluate(VOTE, BOTRYTIS, _profiles(300.0, 100.0))
    assert d.decision == pc.DECISION_UNDETERMINED and d.idiomorph is None
    assert d.note["profile"]["call"] == "MAT1-1"
    assert "disagreement" in d.note["reason"]


def test_inconclusive_profiles_are_undetermined_not_a_flip():
    d = pc.evaluate(VOTE, BOTRYTIS, _profiles(200.0, 210.0))  # 10 bits < 25
    assert d.decision == pc.DECISION_UNDETERMINED
    assert d.note["profile"]["call"] == "inconclusive"


def test_profiles_unavailable_is_undetermined():
    for provider in (None, lambda genes: None, lambda genes: {}):
        d = pc.evaluate(VOTE, BOTRYTIS, provider)
        assert d.decision == pc.DECISION_UNDETERMINED
        assert d.note["profile"] is None and "unavailable" in d.note["reason"]


def test_truncation_at_a_contig_edge_keeps_the_vote():
    stats = [_stat("MAT1-1-1", "MAT1-1", 60.6, 416.0, edge=120), BOTRYTIS[1]]
    d = pc.evaluate(VOTE, stats, _profiles(125.0, 221.0))
    assert d.decision == pc.DECISION_KEPT_EDGE and d.idiomorph is None
    assert d.note["winner_edge_distance"] == 120


def test_two_full_length_genes_do_not_trigger():
    stats = [_stat("MAT1-1-1", "MAT1-1", 99.0, 416.0), _stat("MAT1-2-1", "MAT1-2", 100.0, 326.0)]
    assert pc.evaluate(VOTE, stats, _profiles(125.0, 221.0)) is None


def test_two_partial_genes_do_not_trigger():
    stats = [_stat("MAT1-1-1", "MAT1-1", 60.0, 416.0), _stat("MAT1-2-1", "MAT1-2", 70.0, 326.0)]
    assert pc.evaluate(VOTE, stats, _profiles(125.0, 221.0)) is None


def test_other_side_between_partial_and_full_does_not_trigger():
    stats = [_stat("MAT1-1-1", "MAT1-1", 60.0, 416.0), _stat("MAT1-2-1", "MAT1-2", 85.0, 326.0)]
    assert pc.evaluate(VOTE, stats, _profiles(125.0, 221.0)) is None


def test_unknown_coverage_never_triggers():
    stats = [_stat("MAT1-1-1", "MAT1-1", None, 416.0), _stat("MAT1-2-1", "MAT1-2", 100.0, 326.0)]
    assert pc.evaluate(VOTE, stats, _profiles(125.0, 221.0)) is None


def test_a_winner_with_another_full_core_gene_does_not_trigger():
    stats = BOTRYTIS + [_stat("MAT1-1-9", "MAT1-1", 100.0, 200.0)]
    assert pc.evaluate(VOTE, stats, _profiles(125.0, 221.0)) is None


def test_tie_and_single_candidate_do_not_trigger():
    tie = [{"idiomorph": "MAT1-1", "score": 5.0}, {"idiomorph": "MAT1-2", "score": 5.0}]
    assert pc.evaluate(tie, BOTRYTIS, _profiles(1, 2)) is None
    assert pc.evaluate(VOTE[:1], BOTRYTIS, _profiles(1, 2)) is None


def test_reference_coverage_is_a_capped_length_ratio():
    model = PolishModel("g", None, "core_MAT", "c", 1, 642, "+", [ExonSpan(1, 642)], 90.0, "r", "exonerate_refine")
    assert model_reference_coverage(model, 353) == 60.6
    assert model_reference_coverage(model, 100) == 100.0
    assert model_reference_coverage(model, None) is None


# ---- through run_pipeline ---------------------------------------------------

ORDER = (
    "phylum: P\nloci:\n  - locus_name: MAT\n    vocabulary_type: enum\n"
    "    idiomorph_values: [MAT1-1, MAT1-2]\n    taxonomic_scope: [1]\n    genes:\n"
    "      - {name: MAT1-1-1, role: core_MAT, present_in_idiomorphs: [\"MAT1-1\"]}\n"
    "      - {name: MAT1-2-1, role: core_MAT, present_in_idiomorphs: [\"MAT1-2\"]}\n"
)
KEY = pipeline_module.FamilyKey("P", "MAT")


def _write_records(tmp_path, extra_genes=0):
    (tmp_path / "P").mkdir(exist_ok=True)
    # extra roster genes that are never found push the locus under the fraction floor (relaxed pass)
    (tmp_path / "P" / "order.yml").write_text(ORDER + "".join(
        f"      - {{name: X{i}, role: flanking_conserved}}\n" for i in range(extra_genes)))
    for rid, gene, idiom, aa in (("rec11", "MAT1-1-1", "MAT1-1", 353), ("rec12", "MAT1-2-1", "MAT1-2", 376)):
        d = tmp_path / "P" / "Fam" / rid
        d.mkdir(parents=True, exist_ok=True)
        (d / "metadata.yaml").write_text(
            f"record_id: {rid}\nmating_type: {{locus_name: MAT, idiomorphs: [{idiom}]}}\n")
        (d / "proteins.faa").write_text(f">{rid}|gene_index=0|name={gene}|role=core_MAT\n" + "M" * aa + "\n")


def _hit(gene, rid, start, end, bits):
    return SearchHit(KEY, gene, "core_MAT", "c1", start, end, "+", 90.0, rid, "tblastn_genome", bitscore=bits)


def _model(gene, rid, start, aa, method):
    return PolishModel(gene, KEY, "core_MAT", "c1", start, start + aa * 3 - 1, "+",
                       [ExonSpan(start, start + aa * 3 - 1)], 95.0, rid, method)


def _run(tmp_path, monkeypatch, *, policy, a_start=3000, a_aa=214, a_gene_bits=416.0,
         h_aa=376, profiles=(125.0, 221.0), contig_len=20000, extra_genes=0):
    _write_records(tmp_path, extra_genes)
    (tmp_path / "genome.fa").write_text(">c1\n" + "A" * contig_len + "\n")
    monkeypatch.setattr(pipeline_module, "_translate_model", lambda *a, **k: "MKV")
    monkeypatch.setattr(pc, "load_core_profiles", lambda *a, **k: object() if profiles else None)
    monkeypatch.setattr(pipeline_module, "own_profile_scores",
                        lambda clf, gi, prots: {"MAT1-1-1": profiles[0], "MAT1-2-1": profiles[1]})
    h_start = a_start + a_aa * 3 + 2000
    models = {"MAT1-1-1": (a_start, a_aa, "rec11"), "MAT1-2-1": (h_start, h_aa, "rec12")}

    def polish(method):
        def fn(*, gene_name, window, **kw):
            start, aa, rid = models[gene_name]
            return _model(gene_name, rid, start, aa, method)
        return fn

    hits = [_hit("MAT1-1-1", "rec11", a_start, a_start + a_aa * 3 - 1, a_gene_bits),
            _hit("MAT1-2-1", "rec12", h_start, h_start + h_aa * 3 - 1, 326.0)]
    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None, db_root=tmp_path,
        reference_fasta=tmp_path / "reference.faa", search_localize=lambda *a, **k: hits,
        polish_with_exonerate=polish("exonerate_refine"), polish_with_miniprot=polish("miniprot_refine"),
        partial_core_policy=policy,
    )
    assert len(outcome.results) == 1
    return outcome.results[0]


def test_default_off_keeps_the_bitscore_vote_but_reports_coverage(tmp_path, monkeypatch):
    r = _run(tmp_path, monkeypatch, policy="off")
    assert r.idiomorph == "MAT1-1" and r.idiomorph_vote_note is None
    cov = {e.gene_name: e.coverage for e in r.gene_evidence}
    assert cov == {"MAT1-1-1": 60.6, "MAT1-2-1": 100.0}


def test_default_argument_is_off(tmp_path, monkeypatch):
    import inspect
    assert inspect.signature(run_pipeline).parameters["partial_core_policy"].default == "off"


def test_on_flips_botrytis_when_profiles_agree(tmp_path, monkeypatch):
    r = _run(tmp_path, monkeypatch, policy="coverage_profile")
    assert r.idiomorph == "MAT1-2"
    assert r.idiomorph_vote_note["decision"] == "flipped"
    assert r.confidence != "high" and r.idiomorph_margin is None
    assert [c["idiomorph"] for c in r.idiomorph_candidates][0] == "MAT1-1"  # the vote itself is preserved


def test_on_undetermined_when_profiles_disagree(tmp_path, monkeypatch):
    r = _run(tmp_path, monkeypatch, policy="coverage_profile", profiles=(300.0, 100.0))
    assert r.idiomorph == "undetermined"
    assert r.idiomorph_vote_note["decision"] == "undetermined"
    assert r.idiomorph_vote_note["profile"]["scores"] == {"MAT1-1-1": 300.0, "MAT1-2-1": 100.0}


def test_on_undetermined_when_profiles_missing(tmp_path, monkeypatch):
    r = _run(tmp_path, monkeypatch, policy="coverage_profile", profiles=None)
    assert r.idiomorph == "undetermined" and "unavailable" in r.idiomorph_vote_note["reason"]


def test_truncated_true_mat11_at_the_contig_edge_does_not_flip(tmp_path, monkeypatch):
    r = _run(tmp_path, monkeypatch, policy="coverage_profile", a_start=150)
    assert r.idiomorph == "MAT1-1"
    assert r.idiomorph_vote_note["decision"] == "kept_edge_artefact"


def test_two_full_length_genes_unchanged_under_on(tmp_path, monkeypatch):
    off = _run(tmp_path, monkeypatch, policy="off", a_aa=353)
    on = _run(tmp_path, monkeypatch, policy="coverage_profile", a_aa=353)
    assert on.idiomorph == off.idiomorph == "MAT1-1"
    assert on.idiomorph_vote_note is None and on.confidence == off.confidence
    assert on.idiomorph_margin == off.idiomorph_margin


def test_report_carries_the_vote_note_only_when_set(tmp_path, monkeypatch):
    from MATPredict.detect.report import _result_doc
    off = _result_doc(_run(tmp_path / "a" if (tmp_path / "a").mkdir() is None else None, monkeypatch, policy="off"))
    on = _result_doc(_run(tmp_path / "b" if (tmp_path / "b").mkdir() is None else None, monkeypatch,
                          policy="coverage_profile"))
    assert "idiomorph_vote_note" not in off
    assert on["idiomorph_vote_note"]["decision"] == "flipped"
    assert {g["gene"]: g["coverage"] for g in on["gene_evidence"]} == {"MAT1-1-1": 60.6, "MAT1-2-1": 100.0}


def test_the_relaxed_pass_applies_the_rule_too(tmp_path, monkeypatch):
    off = _run(tmp_path / "a" if (tmp_path / "a").mkdir() is None else None, monkeypatch, policy="off", extra_genes=6)
    on = _run(tmp_path / "b" if (tmp_path / "b").mkdir() is None else None, monkeypatch,
              policy="coverage_profile", extra_genes=6)
    assert off.detection_pass == on.detection_pass == "relaxed"
    assert off.idiomorph == "MAT1-1" and on.idiomorph == "MAT1-2"
    assert on.idiomorph_vote_note["decision"] == "flipped" and off.idiomorph_vote_note is None
