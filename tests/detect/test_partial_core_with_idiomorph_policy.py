"""`--partial-core-policy` and `--idiomorph-policy` are independent and work together.

`--idiomorph-policy` (branch optional-gene-policy) is not on this branch yet, so
these tests skip until it lands; they were run on a scratch merge of both.

Pattern: the Botrytis partial MAT1-1-1 remnant (bitscore 416) beside a full-length
MAT1-2-1 (326), plus an OPTIONAL MAT1-1-3 hit overlapping MAT1-2-1 with the higher
identity (the Parastagonospora failure `core_first` exists for). Expected calls:

    idiomorph-policy  partial-core-policy  call
    legacy            off                  MAT1-1   (MAT1-1-3 supersedes MAT1-2-1; remnant votes)
    legacy            coverage_profile     MAT1-1   (MAT1-2-1 is superseded, so no full-length
                                                     gene is left for this rule: not its business)
    core_first        off                  MAT1-1   (MAT1-2-1 live; the remnant outscores it)
    core_first        coverage_profile     MAT1-2   (the rule's flip)
"""
from __future__ import annotations

import inspect

import pytest

from MATPredict.detect import partial_core as pc
from MATPredict.detect import pipeline as pipeline_module
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.polish import ExonSpan, PolishModel
from MATPredict.detect.search import SearchHit

pytestmark = pytest.mark.skipif(
    "idiomorph_policy" not in inspect.signature(run_pipeline).parameters,
    reason="--idiomorph-policy (origin/optional-gene-policy) is not on this branch")

ORDER = (
    "phylum: P\nloci:\n  - locus_name: MAT\n    vocabulary_type: enum\n"
    "    idiomorph_values: [MAT1-1, MAT1-2]\n    taxonomic_scope: [1]\n    genes:\n"
    "      - {name: MAT1-1-1, role: core_MAT, present_in_idiomorphs: [\"MAT1-1\"]}\n"
    "      - {name: MAT1-2-1, role: core_MAT, present_in_idiomorphs: [\"MAT1-2\"]}\n"
    "      - {name: MAT1-1-3, role: core_MAT, present_in_idiomorphs: [\"MAT1-1\"], optional: true}\n"
)
KEY = pipeline_module.FamilyKey("P", "MAT")


def _records(tmp_path):
    (tmp_path / "P").mkdir(exist_ok=True)
    (tmp_path / "P" / "order.yml").write_text(ORDER)
    for rid, gene, idiom, aa in (("rec11", "MAT1-1-1", "MAT1-1", 353), ("rec12", "MAT1-2-1", "MAT1-2", 376),
                                 ("rec13", "MAT1-1-3", "MAT1-1", 150)):
        d = tmp_path / "P" / "Fam" / rid
        d.mkdir(parents=True, exist_ok=True)
        (d / "metadata.yaml").write_text(f"record_id: {rid}\nmating_type: {{locus_name: MAT, idiomorphs: [{idiom}]}}\n")
        (d / "proteins.faa").write_text(f">{rid}|gene_index=0|name={gene}|role=core_MAT\n" + "M" * aa + "\n")


def _hit(gene, rid, start, end, ident, bits):
    return SearchHit(KEY, gene, "core_MAT", "c1", start, end, "+", ident, rid, "tblastn_genome", bitscore=bits)


def _call(tmp_path, monkeypatch, idiomorph_policy, partial_policy):
    _records(tmp_path)
    (tmp_path / "genome.fa").write_text(">c1\n" + "A" * 20000 + "\n")
    monkeypatch.setattr(pipeline_module, "_translate_model", lambda *a, **k: "MKV")
    monkeypatch.setattr(pc, "load_core_profiles", lambda *a, **k: object())
    monkeypatch.setattr(pipeline_module, "own_profile_scores",
                        lambda clf, gi, prots: {"MAT1-1-1": 125.0, "MAT1-2-1": 221.0})
    h_start = 3000 + 214 * 3 + 2000
    h_end = h_start + 376 * 3 - 1
    models = {"MAT1-1-1": (3000, 214, "rec11"), "MAT1-2-1": (h_start, 376, "rec12")}

    def polish(method):
        def fn(*, gene_name, window, **kw):
            if gene_name not in models:
                return None
            start, aa, rid = models[gene_name]
            return PolishModel(gene_name, KEY, "core_MAT", "c1", start, start + aa * 3 - 1, "+",
                               [ExonSpan(start, start + aa * 3 - 1)], 95.0, rid, method)
        return fn

    hits = [_hit("MAT1-1-1", "rec11", 3000, 3641, 92.7, 416.0),
            _hit("MAT1-2-1", "rec12", h_start, h_end, 99.7, 326.0),
            _hit("MAT1-1-3", "rec13", h_start + 10, h_end - 10, 99.9, 60.0)]  # overlaps MAT1-2-1, higher identity
    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None, db_root=tmp_path,
        reference_fasta=tmp_path / "reference.faa", search_localize=lambda *a, **k: hits,
        polish_with_exonerate=polish("exonerate_refine"), polish_with_miniprot=polish("miniprot_refine"),
        idiomorph_policy=idiomorph_policy, partial_core_policy=partial_policy)
    assert len(outcome.results) == 1
    return outcome.results[0].idiomorph


@pytest.mark.parametrize("idiomorph_policy,partial_policy,expected", [
    ("legacy", "off", "MAT1-1"),
    ("legacy", "coverage_profile", "MAT1-1"),
    ("core_first", "off", "MAT1-1"),
    ("core_first", "coverage_profile", "MAT1-2"),
])
def test_the_four_combinations(tmp_path, monkeypatch, idiomorph_policy, partial_policy, expected):
    assert _call(tmp_path, monkeypatch, idiomorph_policy, partial_policy) == expected
