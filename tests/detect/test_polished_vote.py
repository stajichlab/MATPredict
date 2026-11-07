"""The polished-vote rule (`polished_vote`, `--polished-vote-policy`, #275).

Pure-function tests of `polished_vote.decide`, the same scenarios through `run_pipeline` with stubbed polishing (the scope is pointed at the
test family), then the CLI default and passthrough. The case behind it: Trichoderma viride, where an unpolished 47-bit MAT1-2-1 window beat a
polished 39.7-bit MAT1-1-3 model (truth MAT1-1).
"""
from __future__ import annotations

import argparse
import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

from MATPredict.detect import pipeline as pipeline_module
from MATPredict.detect import polished_vote as pv
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.polish import ExonSpan, PolishModel
from MATPredict.detect.search import SearchHit

KEY = ("Ascomycota", "MAT")
FAM = Family(
    FamilyKey("Ascomycota", "MAT"), "enum", ["MAT1-1", "MAT1-2"], None,
    [
        {"name": "SLA2", "role": "flanking_conserved"},
        {"name": "MAT1-1-1", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-1"]},
        {"name": "MAT1-1-3", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-1"]},
        {"name": "MAT1-2-1", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-2"]},
    ],
    [4890],
)


def G(gene, status, bits):
    return pv.GeneStat(gene, status, bits)


def test_default_is_off_and_check_policy():
    assert pv.DEFAULT_POLISHED_VOTE_POLICY == "off"
    assert inspect.signature(run_pipeline).parameters["polished_vote_policy"].default == "off"
    for p in pv.POLISHED_VOTE_POLICIES:
        assert pv.check_policy(p) == p
    with pytest.raises(ValueError):
        pv.check_policy("bogus")


def test_off_never_flips():
    stats = [G("MAT1-2-1", "unpolished", 47.0), G("MAT1-1-3", "polished_single", 39.7)]
    assert pv.decide(FAM, KEY, "MAT1-2", stats, policy="off") is None


def test_a_polished_idiomorph_beats_an_unpolished_one_the_viride_case():
    stats = [G("MAT1-2-1", "unpolished", 47.0), G("MAT1-1-3", "polished_single", 39.7)]
    d = pv.decide(FAM, KEY, "MAT1-2", stats)
    assert d.idiomorph == "MAT1-1"
    assert d.note["rule"] == "polished_vote" and d.note["decision"] == "flipped" and d.note["from"] == "MAT1-2"
    assert d.note["set_aside_gene"] == "MAT1-2-1" and d.note["polished_gene"] == "MAT1-1-3" and d.note["polished_bitscore"] == 39.7


def test_nothing_changes_when_the_called_idiomorph_has_a_polished_gene():
    stats = [G("MAT1-2-1", "polished_disagree", 47.0), G("MAT1-1-3", "polished_single", 39.7)]
    assert pv.decide(FAM, KEY, "MAT1-2", stats) is None


def test_nothing_changes_when_the_other_idiomorph_is_unpolished_too():
    stats = [G("MAT1-2-1", "unpolished", 47.0), G("MAT1-1-3", "unpolished", 60.0)]
    assert pv.decide(FAM, KEY, "MAT1-2", stats) is None


def test_nothing_changes_when_only_one_idiomorph_has_evidence():
    assert pv.decide(FAM, KEY, "MAT1-2", [G("MAT1-2-1", "unpolished", 47.0)]) is None


def test_the_best_polished_gene_of_the_other_idiomorph_is_named():
    stats = [G("MAT1-2-1", "unpolished", 47.0), G("MAT1-1-3", "polished_single", 39.7), G("MAT1-1-1", "polished_agree", 120.0)]
    d = pv.decide(FAM, KEY, "MAT1-2", stats)
    assert d.note["polished_gene"] == "MAT1-1-1"


def test_flank_genes_and_lost_genes_do_not_vote_and_scope_is_ascomycota_mat_only():
    stats = [G("MAT1-2-1", "unpolished", 47.0), G("SLA2", "polished_agree", 1700.0)]
    assert pv.decide(FAM, KEY, "MAT1-2", stats) is None                  # SLA2 declares no idiomorph
    stats = [G("MAT1-2-1", "unpolished", 47.0), G("MAT1-1-3", "polished_single", 39.7)]
    assert pv.decide(FAM, ("Ascomycota", "MATsc"), "MAT1-2", stats) is None
    assert pv.decide(FAM, ("Basidiomycota", "MAT"), "MAT1-2", stats) is None
    assert pv.decide(FAM, KEY, "undetermined", stats) is None


# ---- through run_pipeline ---------------------------------------------------------------------------------------------------------

ORDER = (
    "phylum: P\nloci:\n  - locus_name: MATX\n    vocabulary_type: enum\n"
    "    idiomorph_values: [alpha, a]\n    taxonomic_scope: [1]\n    genes:\n"
    "      - {name: ALPHA1, role: core_MAT, present_in_idiomorphs: [\"alpha\"]}\n"
    "      - {name: A1, role: core_MAT, present_in_idiomorphs: [\"a\"]}\n"
    "      - {name: FL1, role: flanking_conserved}\n"
)
PKEY = pipeline_module.FamilyKey("P", "MATX")


def _write_records(tmp_path):
    (tmp_path / "P").mkdir(exist_ok=True)
    (tmp_path / "P" / "order.yml").write_text(ORDER)
    for rid, gene, idiom, role in (("r1", "ALPHA1", "alpha", "core_MAT"), ("r3", "A1", "a", "core_MAT"), ("r4", "FL1", "alpha", "flanking_conserved")):
        d = tmp_path / "P" / "Fam" / rid
        d.mkdir(parents=True, exist_ok=True)
        (d / "metadata.yaml").write_text(f"record_id: {rid}\nmating_type: {{locus_name: MATX, idiomorphs: [{idiom}]}}\n")
        (d / "proteins.faa").write_text(f">{rid}|gene_index=0|name={gene}|role={role}\n" + "M" * 300 + "\n")


def _run(tmp_path, monkeypatch, *, policy):
    """ALPHA1: a polished model scoring 39.7; A1: an unpolished tblastn hit scoring 47 (the opposite idiomorph)."""
    monkeypatch.setattr(pv, "SCOPE_FAMILY", ("P", "MATX"))
    _write_records(tmp_path)
    (tmp_path / "genome.fa").write_text(">c1\n" + "A" * 60000 + "\n")
    monkeypatch.setattr(pipeline_module, "_translate_model", lambda *a, **k: "MKV")
    aa = 290
    a_start, b_start = 30000, 33000
    hits = [
        SearchHit(PKEY, "ALPHA1", "core_MAT", "c1", a_start, a_start + aa * 3 - 1, "+", 40.0, "r1", "tblastn_genome", bitscore=39.7),
        SearchHit(PKEY, "A1", "core_MAT", "c1", b_start, b_start + 100 * 3 - 1, "+", 45.0, "r3", "tblastn_genome", bitscore=47.0),
        SearchHit(PKEY, "FL1", "flanking_conserved", "c1", 20000, 20000 + aa * 3 - 1, "+", 90.0, "r4", "tblastn_genome", bitscore=300.0),
    ]

    def polish(*, gene_name, window, **kw):
        if gene_name == "ALPHA1":
            return PolishModel("ALPHA1", PKEY, "core_MAT", "c1", a_start, a_start + aa * 3 - 1, "+",
                               [ExonSpan(a_start, a_start + aa * 3 - 1)], 40.0, "r1", "miniprot_refine")
        if gene_name == "FL1":
            return PolishModel("FL1", PKEY, "flanking_conserved", "c1", 20000, 20000 + aa * 3 - 1, "+",
                               [ExonSpan(20000, 20000 + aa * 3 - 1)], 90.0, "r4", "miniprot_refine")
        return None

    out = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None, db_root=tmp_path,
        reference_fasta=tmp_path / "reference.faa", search_localize=lambda *a, **k: hits,
        polish_with_exonerate=polish, polish_with_miniprot=polish, polished_vote_policy=policy,
    )
    return out.results


def test_the_vote_through_run_pipeline(tmp_path, monkeypatch):
    (tmp_path / "a").mkdir(); (tmp_path / "b").mkdir()
    off = _run(tmp_path / "a", monkeypatch, policy="off")
    on = _run(tmp_path / "b", monkeypatch, policy="on")
    assert len(off) == 1 and len(on) == 1, "the stub locus must be admitted (a polished flank gene is needed)"
    assert off[0].idiomorph == "a"                     # the unpolished 47-bit fragment wins the bitscore vote
    assert on[0].idiomorph == "alpha"                  # the polished model wins under the rule
    assert on[0].idiomorph_vote_note["rule"] == "polished_vote" and on[0].idiomorph_vote_note["from"] == "a"
    assert on[0].idiomorph_margin is None and on[0].confidence != "high"
    assert off[0].idiomorph_vote_note is None


# ---- CLI ---------------------------------------------------------------------------------------------------------------------------

def test_cli_default_and_passthrough(monkeypatch, tmp_path):
    from MATPredict.detect import cli as detect_cli
    from MATPredict.detect.pipeline import DetectionOutcome

    assert detect_cli._polished_vote_policy_from_args(argparse.Namespace()) == "off"
    monkeypatch.setattr(detect_cli, "MatpredictConfig", SimpleNamespace(from_env=lambda repo_root: SimpleNamespace(db_root=Path("db"))))
    seen: list = []

    def fake_run_pipeline(**kwargs):
        seen.append(kwargs.get("polished_vote_policy", "MISSING"))
        return DetectionOutcome(results=[])

    monkeypatch.setattr(detect_cli, "build_reference_fasta", lambda db_root, out, family_keys=None, **kw: out)
    monkeypatch.setattr(detect_cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(detect_cli, "write_detection_gff3", lambda *a, **k: None)
    monkeypatch.setattr(detect_cli, "write_detection_report", lambda outcome, path, run=None: None)
    base = dict(genome="g.fa", proteins=None, taxid=None, out_dir=str(tmp_path / "out"), evidence_diagnostics=None, min_hits=1,
                min_identity=None, require_core_role=False, exclude_records="", emit_cds_fasta=False, phylum="Mucoromycota")
    assert detect_cli._cmd_detect(SimpleNamespace(**base)) == 0
    assert detect_cli._cmd_detect(SimpleNamespace(**base, polished_vote_policy="on")) == 0
    assert seen == ["off", "on"]


def test_no_preset_turns_the_vote_rule_on():
    from MATPredict.detect.cli import PRESETS
    assert all("polished_vote_policy" not in v for v in PRESETS.values())
