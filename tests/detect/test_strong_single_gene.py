"""`--strong-single-gene-policy`: a lone, strong core gene is a low-confidence locus.

Curator's decisions 2026-10-17. Pattern (three Parastagonospora genomes): the true
MAT1-2-1 is a lone modelled gene at 99%, no flank, no second modelled gene; the
evidence floor, the two-modelled-gene bar and the fraction floor each drop it.
With the rule ON (default OFF) it is reported as `idiomorph_gene_only` at `low`;
an HMG paralog at 56% is not; two qualifying loci of opposite idiomorph are both
reported and flagged through `two_idiomorphs`.
"""
from __future__ import annotations

import argparse
import inspect

import pytest
import yaml

from MATPredict.detect import pipeline as pipeline_module
from MATPredict.detect.cli import _strong_single_gene_policy_from_args, register_subcommands
from MATPredict.detect.idiomorph import LOCUS_CLASS_IDIOMORPH_ONLY
from MATPredict.detect.pipeline import DetectionResult, GeneEvidence, run_pipeline
from MATPredict.detect.polish import ExonSpan, PolishModel
from MATPredict.detect.report import write_detection_report
from MATPredict.detect.search import SearchHit
from MATPredict.detect.strong_single_gene import (
    DEFAULT_STRONG_SINGLE_GENE_POLICY, StrongSingleGeneParams, evaluate_strong_single_gene,
)

ORDER = (
    "phylum: P\nloci:\n  - locus_name: MAT\n    vocabulary_type: enum\n"
    "    idiomorph_values: [MAT1-1, MAT1-2]\n    taxonomic_scope: [1]\n    genes:\n"
    "      - {name: MAT1-1-1, role: core_MAT, present_in_idiomorphs: [\"MAT1-1\"]}\n"
    "      - {name: MAT1-2-1, role: core_MAT, present_in_idiomorphs: [\"MAT1-2\"]}\n"
    "      - {name: APN2, role: flanking_conserved}\n"
    "      - {name: SLA2, role: flanking_conserved}\n"
)
KEY = pipeline_module.FamilyKey("P", "MAT")
REF_AA = 300


def _records(tmp_path):
    (tmp_path / "P").mkdir(exist_ok=True)
    (tmp_path / "P" / "order.yml").write_text(ORDER)
    for rid, gene, idiom in (("rec11", "MAT1-1-1", "MAT1-1"), ("rec12", "MAT1-2-1", "MAT1-2")):
        d = tmp_path / "P" / "Fam" / rid
        d.mkdir(parents=True, exist_ok=True)
        (d / "metadata.yaml").write_text(
            f"record_id: {rid}\nmating_type: {{locus_name: MAT, idiomorphs: [{idiom}]}}\n")
        (d / "proteins.faa").write_text(
            f">{rid}|gene_index=0|name={gene}|role=core_MAT\n" + "M" * REF_AA + "\n")
    (tmp_path / "genome.fa").write_text(
        ">c1\n" + "A" * 40000 + "\n>c2\n" + "A" * 40000 + "\n")


def _hit(gene, contig, start, ident, role="core_MAT", aa=REF_AA, bits=300.0):
    rid = {"MAT1-1-1": "rec11"}.get(gene, "rec12")
    return SearchHit(KEY, gene, role, contig, start, start + aa * 3 - 1, "+", ident, rid,
                     "tblastn_genome", bitscore=bits)


def _run(tmp_path, hits, *, models, policy, **kw):
    """`models`: gene -> (identity, aa) the polishing tools build for it; absent = no model."""
    _records(tmp_path)

    def polish(method):
        def fn(*, gene_name, window, **k):
            if gene_name not in models:
                return None
            ident, aa = models[gene_name]
            contig, w_start, w_end = window
            h = next(x for x in hits if x.gene_name == gene_name and x.contig == contig
                     and x.start <= w_end and x.end >= w_start)
            rid = "rec11" if gene_name == "MAT1-1-1" else "rec12"
            role = "core_MAT" if gene_name.startswith("MAT") else "flanking_conserved"
            return PolishModel(gene_name, KEY, role, h.contig, h.start, h.start + aa * 3 - 1, "+",
                               [ExonSpan(h.start, h.start + aa * 3 - 1)], ident, rid, method)
        return fn

    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None, db_root=tmp_path,
        reference_fasta=tmp_path / "reference.faa", search_localize=lambda *a, **k: hits,
        polish_with_exonerate=polish("exonerate_refine"), polish_with_miniprot=polish("miniprot_refine"),
        strong_single_gene_policy=policy, **kw)


LONE = [_hit("MAT1-2-1", "c1", 10000, 99.4)]
LONE_MODEL = {"MAT1-2-1": (99.4, REF_AA)}


def test_the_default_is_off_and_nothing_is_reported(tmp_path):
    assert DEFAULT_STRONG_SINGLE_GENE_POLICY == "off"
    outcome = _run(tmp_path, LONE, models=LONE_MODEL, policy="off")
    assert outcome.results == []
    assert [n.family_key for n in outcome.not_detected] == [KEY]


def test_on_reports_the_lone_strong_gene_as_idiomorph_gene_only_at_low(tmp_path):
    outcome = _run(tmp_path, LONE, models=LONE_MODEL, policy="on")
    [r] = outcome.results
    assert r.idiomorph == "MAT1-2"
    assert r.confidence == "low"
    assert r.locus_class == LOCUS_CLASS_IDIOMORPH_ONLY
    assert r.strong_single_gene["core_gene"] == "MAT1-2-1"
    assert r.strong_single_gene["core_identity"] == 99.4
    assert r.strong_single_gene["core_coverage"] == 100.0
    assert r.strong_single_gene["cluster_identity"] == 99.4
    assert outcome.not_detected == []
    assert outcome.two_idiomorphs == []

    path = tmp_path / "report.yaml"
    write_detection_report(outcome, path)
    row = yaml.safe_load(path.read_text())["detected"][0]
    assert (row["confidence"], row["locus_class"]) == ("low", "idiomorph_gene_only")
    assert row["strong_single_gene"]["min_core_identity"] == 85.0


def test_a_superseded_cross_hit_in_the_cluster_does_not_make_it_a_second_gene(tmp_path):
    # the Parastagonospora cluster: a MAT1-1-1 cross-hit over the same gene, weaker
    hits = LONE + [_hit("MAT1-1-1", "c1", 10010, 45.0, aa=290)]
    outcome = _run(tmp_path, hits, models=LONE_MODEL, policy="on")
    [r] = outcome.results
    assert r.idiomorph == "MAT1-2" and r.strong_single_gene is not None


@pytest.mark.parametrize("ident", [56.0, 70.0, 84.9])
def test_a_paralog_like_hit_below_the_gate_is_not_reported(tmp_path, ident):
    outcome = _run(tmp_path, [_hit("MAT1-2-1", "c1", 10000, ident)],
                   models={"MAT1-2-1": (ident, REF_AA)}, policy="on")
    assert outcome.results == []


def test_the_gate_is_inclusive_at_85_and_uses_the_model_identity(tmp_path):
    ok = _run(tmp_path, [_hit("MAT1-2-1", "c1", 10000, 85.0)],
              models={"MAT1-2-1": (85.0, REF_AA)}, policy="on")
    assert len(ok.results) == 1
    # a strong raw hit does not stand for the locus when the MODEL is below the gate
    no = _run(tmp_path, [_hit("MAT1-2-1", "c1", 10000, 99.0)],
              models={"MAT1-2-1": (80.0, REF_AA)}, policy="on")
    assert no.results == []


def test_a_short_model_is_not_near_full_length(tmp_path):
    outcome = _run(tmp_path, LONE, models={"MAT1-2-1": (99.4, int(REF_AA * 0.6))}, policy="on")
    assert outcome.results == []


def test_an_unmodelled_hit_is_not_reported(tmp_path):
    outcome = _run(tmp_path, LONE, models={}, policy="on")
    assert outcome.results == []


def test_a_flanking_gene_alone_is_not_reported(tmp_path):
    hits = [_hit("APN2", "c1", 10000, 99.0, role="flanking_conserved")]
    outcome = _run(tmp_path, hits, models={}, policy="on")
    assert outcome.results == []


def test_opposite_idiomorph_loci_are_both_reported_and_flagged(tmp_path):
    hits = [_hit("MAT1-1-1", "c1", 10000, 94.5), _hit("MAT1-2-1", "c2", 12000, 99.0)]
    models = {"MAT1-1-1": (94.5, REF_AA), "MAT1-2-1": (99.0, REF_AA)}
    outcome = _run(tmp_path, hits, models=models, policy="on")
    assert sorted(r.idiomorph for r in outcome.results) == ["MAT1-1", "MAT1-2"]
    assert all(r.confidence == "low" and r.locus_class == LOCUS_CLASS_IDIOMORPH_ONLY
               for r in outcome.results)
    [st] = outcome.two_idiomorphs
    assert st["arrangement"] == "unlinked"
    assert "homothallism" not in st["supported_causes"]
    assert {st["plus_call"]["idiomorph"], st["minus_call"]["idiomorph"]} == {"MAT1-1", "MAT1-2"}
    # off: neither is reported
    off = _run(tmp_path, hits, models=models, policy="off")
    assert off.results == [] and off.two_idiomorphs == []


def test_two_copies_of_one_idiomorph_report_the_best_only(tmp_path):
    hits = [_hit("MAT1-2-1", "c1", 10000, 90.0), _hit("MAT1-2-1", "c2", 12000, 99.0)]
    outcome = _run(tmp_path, hits, models={"MAT1-2-1": (99.0, REF_AA)}, policy="on")
    assert len(outcome.results) == 1
    assert len({r.idiomorph for r in outcome.results}) == 1
    assert outcome.two_idiomorphs == []


def test_a_family_that_is_already_called_is_left_alone(tmp_path):
    # core + two flanks on c1 (a normal call) and a lone strong MAT1-1-1 on c2
    hits = [_hit("MAT1-2-1", "c1", 10000, 99.0), _hit("APN2", "c1", 12000, 90.0, "flanking_conserved", aa=100),
            _hit("SLA2", "c1", 14000, 90.0, "flanking_conserved", aa=100),
            _hit("MAT1-1-1", "c2", 12000, 99.0)]
    models = {"MAT1-2-1": (99.0, REF_AA), "MAT1-1-1": (99.0, REF_AA),
              "APN2": (90.0, 100), "SLA2": (90.0, 100)}
    off = _run(tmp_path, hits, models=models, policy="off")
    on = _run(tmp_path, hits, models=models, policy="on")
    sig = lambda o: [(r.contig, r.start, r.end, r.idiomorph, r.confidence, r.locus_class) for r in o.results]
    assert len(off.results) == 1 and off.results[0].contig == "c1"
    assert sig(on) == sig(off)
    assert all(r.strong_single_gene is None for r in on.results)


def test_it_coexists_with_the_partial_core_policy(tmp_path):
    from MATPredict.detect import partial_core as pc
    outcome = _run(tmp_path, LONE, models=LONE_MODEL, policy="on",
                   partial_core_policy=pc.POLICY_COVERAGE_PROFILE)
    [r] = outcome.results
    assert r.idiomorph == "MAT1-2" and r.confidence == "low"


@pytest.mark.skipif("idiomorph_policy" not in inspect.signature(run_pipeline).parameters,
                    reason="--idiomorph-policy (origin/optional-gene-policy, PR #89) is not on this branch")
@pytest.mark.parametrize("policy", ["legacy", "core_first"])
def test_it_coexists_with_the_idiomorph_policy(tmp_path, policy):
    # the Parastagonospora cluster: an optional-style cross-hit overlapping the true gene
    hits = LONE + [_hit("MAT1-1-1", "c1", 10010, 45.0, aa=290)]
    outcome = _run(tmp_path, hits, models=LONE_MODEL, policy="on", idiomorph_policy=policy)
    [r] = outcome.results
    assert r.idiomorph == "MAT1-2" and r.strong_single_gene is not None


def test_the_policy_is_a_cli_flag_default_off(tmp_path):
    parser = argparse.ArgumentParser()
    register_subcommands(parser.add_subparsers(dest="cmd"))
    base = ["detect", "--genome", "g.fa", "--out", "o"]
    assert parser.parse_args(base).strong_single_gene_policy == "off"
    assert parser.parse_args(base + ["--strong-single-gene-policy", "on"]).strong_single_gene_policy == "on"
    assert _strong_single_gene_policy_from_args(argparse.Namespace()) == "off"


def _candidate(identity=99.4, coverage=100.0, role="core_MAT", status="polished_agree",
               idiomorph="MAT1-2", extra=()):
    ev = [GeneEvidence("MAT1-2-1", role, "c1", 1, 900, "+", identity, coverage, "rec12",
                       "miniprot_refine", status=status), *extra]
    return DetectionResult(family_key=KEY, contig="c1", start=1, end=900, confidence="low",
                           idiomorph=idiomorph, ambiguous_with=[], genes_found=["MAT1-2-1"],
                           genes_missing=[], fragmented=False, gene_evidence=ev, polished_genes=1)


def test_evaluate_unit():
    assert evaluate_strong_single_gene(_candidate()) is not None
    assert evaluate_strong_single_gene(_candidate(identity=84.9)) is None
    assert evaluate_strong_single_gene(_candidate(coverage=89.9)) is None
    assert evaluate_strong_single_gene(_candidate(coverage=None)) is None
    assert evaluate_strong_single_gene(_candidate(role="flanking_conserved")) is None
    assert evaluate_strong_single_gene(_candidate(status="unpolished")) is None
    assert evaluate_strong_single_gene(_candidate(idiomorph="undetermined")) is None
    second = GeneEvidence("APN2", "flanking_conserved", "c1", 2000, 2500, "+", 70.0, 90.0, "rec12",
                          "tblastn_genome")
    assert evaluate_strong_single_gene(_candidate(extra=[second])) is None
    lost = GeneEvidence("MAT1-1-1", "core_MAT", "c1", 5, 800, "+", 45.0, 90.0, "rec11",
                        "tblastn_genome", lost_to="MAT1-2-1")
    assert evaluate_strong_single_gene(_candidate(extra=[lost])) is not None
    assert evaluate_strong_single_gene(
        _candidate(identity=80.0), params=StrongSingleGeneParams(min_core_identity=80.0)) is not None


def test_on_never_changes_a_result_the_rule_did_not_make(tmp_path):
    # core + one flank: admitted by the ordinary floor; the rule must leave it as `off` has it
    hits = [_hit("MAT1-2-1", "c1", 10000, 99.0), _hit("APN2", "c1", 12000, 90.0, "flanking_conserved", aa=100),
            _hit("MAT1-1-1", "c2", 12000, 99.0)]
    models = {"MAT1-2-1": (99.0, REF_AA), "MAT1-1-1": (99.0, REF_AA), "APN2": (90.0, 100)}
    off = _run(tmp_path, hits, models=models, policy="off")
    on = _run(tmp_path, hits, models=models, policy="on")
    sig = lambda o: sorted((r.contig, r.start, r.end, r.idiomorph, r.confidence, r.locus_class) for r in o.results)
    extra = [r for r in on.results if r.strong_single_gene is None]
    assert sorted((r.contig, r.start, r.end, r.idiomorph, r.confidence, r.locus_class) for r in extra) == sig(off)
    assert [(n.family_key, n.reason) for n in on.not_detected] == [(n.family_key, n.reason) for n in off.not_detected]
