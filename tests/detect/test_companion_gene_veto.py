"""Companion-gene veto of the partial-coverage rule (curator 2026-10-18).

When the rule would move the call away from X and X's declared companion gene (family
`idiomorph_companions`, e.g. MAT1-1: [MAT1-1-5]) is intact in the locus (modelled, >= 90% of
its reference), the locus holds both idiomorphs: `undetermined`, reason in `idiomorph_vote_note`.
"""
from __future__ import annotations

import pytest

from MATPredict.detect import partial_core as pc
from MATPredict.detect import pipeline as pipeline_module
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.polish import ExonSpan, PolishModel
from MATPredict.detect.search import SearchHit
from tests.detect.test_partial_core_vote import BOTRYTIS, VOTE, _profiles

COMP = {"MAT1-1": ["MAT1-1-5"]}


def _c(cov, bits=389.0, gene="MAT1-1-5"):
    return [pc.CoreStat(gene, "", cov, bits, None)]


def test_intact_companion_makes_it_undetermined_with_the_reason():
    d = pc.evaluate(VOTE, BOTRYTIS, _profiles(125.0, 221.0), companions=COMP, companion_stats=_c(100.0))
    assert d.decision == pc.DECISION_UNDETERMINED and d.idiomorph is None
    assert d.note["veto"] == "companion_intact" and d.note["companion_gene"] == "MAT1-1-5"
    assert "companion gene MAT1-1-5 intact" in d.note["reason"]


@pytest.mark.parametrize("stats", [[], _c(60.0), _c(89.9), _c(None)])
def test_absent_or_partial_companion_flips_as_before(stats):
    d = pc.evaluate(VOTE, BOTRYTIS, _profiles(125.0, 221.0), companions=COMP, companion_stats=stats)
    assert d.decision == pc.DECISION_FLIPPED and d.idiomorph == "MAT1-2"


def test_no_companions_declared_is_unchanged():
    for comps in (None, {}, {"MAT1-2": ["MAT1-1-5"]}):  # the last declares it for the other idiomorph
        d = pc.evaluate(VOTE, BOTRYTIS, _profiles(125.0, 221.0), companions=comps, companion_stats=_c(100.0))
        assert d.decision == pc.DECISION_FLIPPED


def test_a_companion_of_the_other_idiomorph_does_not_veto():
    d = pc.evaluate(VOTE, BOTRYTIS, _profiles(125.0, 221.0), companions={"MAT1-2": ["MAT1-2-10"]},
                    companion_stats=_c(100.0, gene="MAT1-2-10"))
    assert d.decision == pc.DECISION_FLIPPED


def test_veto_does_not_hide_the_edge_guard_or_a_non_trigger():
    edge = [pc.CoreStat("MAT1-1-1", "MAT1-1", 60.6, 416.0, 100), BOTRYTIS[1]]
    d = pc.evaluate(VOTE, edge, _profiles(125.0, 221.0), companions=COMP, companion_stats=_c(100.0))
    assert d.decision == pc.DECISION_KEPT_EDGE
    full = [pc.CoreStat("MAT1-1-1", "MAT1-1", 99.0, 416.0, 5000), BOTRYTIS[1]]
    assert pc.evaluate(VOTE, full, _profiles(1, 2), companions=COMP, companion_stats=_c(100.0)) is None


# ---- through run_pipeline ---------------------------------------------------

KEY = pipeline_module.FamilyKey("P", "MAT")
ORDER = (
    "phylum: P\nloci:\n  - locus_name: MAT\n    vocabulary_type: enum\n"
    "    idiomorph_values: [MAT1-1, MAT1-2]\n    taxonomic_scope: [1]\n"
    "{companions}    genes:\n"
    "      - {{name: MAT1-1-1, role: core_MAT, present_in_idiomorphs: [\"MAT1-1\"]}}\n"
    "      - {{name: MAT1-2-1, role: core_MAT, present_in_idiomorphs: [\"MAT1-2\"]}}\n"
    "      - {{name: MAT1-1-5, role: core_MAT, present_in_idiomorphs: [\"MAT1-1\"], optional: true}}\n"
)
GENES = (("rec11", "MAT1-1-1", "MAT1-1", 353), ("rec12", "MAT1-2-1", "MAT1-2", 376), ("rec15", "MAT1-1-5", "MAT1-1", 260))


def _run(tmp_path, monkeypatch, *, declared, c_aa=260):
    comp = "    idiomorph_companions: {MAT1-1: [MAT1-1-5]}\n" if declared else ""
    (tmp_path / "P").mkdir()
    (tmp_path / "P" / "order.yml").write_text(ORDER.format(companions=comp))
    for rid, gene, idiom, aa in GENES:
        d = tmp_path / "P" / "Fam" / rid
        d.mkdir(parents=True)
        (d / "metadata.yaml").write_text(f"record_id: {rid}\nmating_type: {{locus_name: MAT, idiomorphs: [{idiom}]}}\n")
        (d / "proteins.faa").write_text(f">{rid}|gene_index=0|name={gene}|role=core_MAT\n" + "M" * aa + "\n")
    (tmp_path / "genome.fa").write_text(">c1\n" + "A" * 30000 + "\n")
    monkeypatch.setattr(pipeline_module, "_translate_model", lambda *a, **k: "MKV")
    monkeypatch.setattr(pc, "load_core_profiles", lambda *a, **k: object())
    monkeypatch.setattr(pipeline_module, "own_profile_scores",
                        lambda clf, gi, prots: {"MAT1-1-1": 125.0, "MAT1-2-1": 221.0})
    models = {"MAT1-1-1": (3000, 214, "rec11"), "MAT1-2-1": (6000, 376, "rec12"), "MAT1-1-5": (9000, c_aa, "rec15")}

    def polish(method):
        def fn(*, gene_name, window, **kw):
            start, aa, rid = models[gene_name]
            return PolishModel(gene_name, KEY, "core_MAT", "c1", start, start + aa * 3 - 1, "+",
                               [ExonSpan(start, start + aa * 3 - 1)], 95.0, rid, method)
        return fn

    hits = [SearchHit(KEY, g, "core_MAT", "c1", s, s + aa * 3 - 1, "+", 90.0, rid, "tblastn_genome", bitscore=b)
            for (g, rid, s, aa, b) in (("MAT1-1-1", "rec11", 3000, 214, 416.0), ("MAT1-2-1", "rec12", 6000, 376, 326.0),
                                       ("MAT1-1-5", "rec15", 9000, c_aa, 389.0))]
    out = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None, db_root=tmp_path,
        reference_fasta=tmp_path / "reference.faa", search_localize=lambda *a, **k: hits,
        polish_with_exonerate=polish("exonerate_refine"), polish_with_miniprot=polish("miniprot_refine"),
        partial_core_policy="coverage_profile")
    assert len(out.results) == 1
    return out.results[0]


def test_roster_key_is_read(tmp_path, monkeypatch):
    from MATPredict.detect.family_registry import load_all_families
    _run(tmp_path, monkeypatch, declared=True)
    assert load_all_families(tmp_path)[0].idiomorph_companions == {"MAT1-1": ["MAT1-1-5"]}


def test_pipeline_companion_intact_is_undetermined(tmp_path, monkeypatch):
    r = _run(tmp_path, monkeypatch, declared=True)
    assert r.idiomorph == "undetermined"
    assert r.idiomorph_vote_note["veto"] == "companion_intact"


def test_pipeline_partial_companion_still_flips(tmp_path, monkeypatch):
    r = _run(tmp_path, monkeypatch, declared=True, c_aa=100)
    assert r.idiomorph == "MAT1-2" and r.idiomorph_vote_note["decision"] == "flipped"


def test_pipeline_nothing_declared_flips_as_before(tmp_path, monkeypatch):
    r = _run(tmp_path, monkeypatch, declared=False)
    assert r.idiomorph == "MAT1-2"


def test_the_ascomycota_roster_declares_only_mat11_5():
    from pathlib import Path
    from MATPredict.detect.family_registry import load_all_families
    fams = {f.key: f for f in load_all_families(Path(__file__).resolve().parents[2] / "db")}
    comp = {k: f.idiomorph_companions for k, f in fams.items() if f.idiomorph_companions}
    assert comp == {pipeline_module.FamilyKey("Ascomycota", "MAT"): {"MAT1-1": ["MAT1-1-5"]}}
