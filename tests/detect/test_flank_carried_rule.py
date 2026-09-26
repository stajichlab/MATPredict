"""A call whose core genes were never modelled rests on its flanks alone.

Curator's ruling 2026-09-26 (docs/notes/2026-09-26_cauris-flank-carried-and-
calbicans-zygosity.md; measured in docs/notes/2026-09-26_polish-cap-measured-
and-serinales-scan.md). In the Serinales-wide scan 131 of 2,647 calls had no
modelled core gene. 52 were Debaryomyces artefacts: the PAP1-OBP1-PIK1 block
sits ~700 kb from the real MTL genes, and a 29% MTLA2 fragment 8 kb away made
a spurious call. The rule:

1. core hit inside the flank span (+-3 kb): keep, cap at `low`, class
   `partial_locus`, flag `idiomorph_unmodelled: true`;
2. core hit outside it: withhold, like any bar failure.

Generic over the family's own flanking genes (role `flanking_*`), not
hardcoded to PAP1/OBP1/PIK1.
"""
import yaml

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.flank_carried import (
    FLANK_SPAN_PADDING_BP, WITHHELD_FLANK_CARRIED, apply_flank_carried_rule,
)
from MATPredict.detect.idiomorph import LOCUS_CLASS_MAT, LOCUS_CLASS_PARTIAL
from MATPredict.detect.pipeline import DetectionResult, GeneEvidence, run_pipeline
from MATPredict.detect.report import write_detection_report
from MATPredict.detect.search import SearchHit

from tests.detect.test_pipeline import _no_polish, _write_order, _write_record

KEY = FamilyKey("P", "aLocus")


def _ev(gene, role, start, end, status, method="tblastn_genome", contig="c1"):
    return GeneEvidence(gene, role, contig, start, end, "+", 40.0, 50.0, "rec1",
                        method, status=status)


def _result(evidence, confidence="medium", locus_class=LOCUS_CLASS_MAT):
    return DetectionResult(
        family_key=KEY, contig="c1", start=min(e.start for e in evidence),
        end=max(e.end for e in evidence), confidence=confidence, idiomorph="a1",
        ambiguous_with=[], genes_found=sorted({e.gene_name for e in evidence}),
        genes_missing=[], fragmented=False, gene_evidence=evidence,
        locus_class=locus_class, polished_genes=2,
    )


FLANKS = [
    _ev("flk1", "flanking_conserved", 1_000, 2_000, "polished_single"),
    _ev("flk2", "flanking_variable", 10_000, 11_000, "polished_agree"),
]


def test_the_padding_is_three_kb():
    assert FLANK_SPAN_PADDING_BP == 3_000


def test_a_call_with_a_modelled_core_gene_is_untouched():
    r = _result(FLANKS + [_ev("mfa1", "core_MAT", 50_000, 51_000, "polished_single")])
    kept, withheld = apply_flank_carried_rule([r])
    assert kept == [r] and withheld == []


def test_an_annotated_core_gene_counts_as_modelled():
    r = _result(FLANKS + [_ev("mfa1", "core_MAT", 50_000, 51_000, "not_polish_candidate",
                              method="diamond_proteome")])
    kept, withheld = apply_flank_carried_rule([r])
    assert kept == [r] and withheld == []


def test_an_unmodelled_core_inside_the_flank_span_is_kept_at_low():
    r = _result(FLANKS + [_ev("mfa1", "core_MAT", 5_000, 5_500, "unpolished")],
                confidence="high")
    [kept], withheld = apply_flank_carried_rule([r])
    assert withheld == []
    assert kept.confidence == "low"
    assert kept.locus_class == LOCUS_CLASS_PARTIAL
    assert kept.idiomorph_unmodelled is True


def test_the_padding_extends_the_span_on_both_sides():
    just_inside = _result(FLANKS + [_ev("mfa1", "core_MAT", 13_500, 14_000, "unpolished")])
    just_outside = _result(FLANKS + [_ev("mfa1", "core_MAT", 13_500, 14_001, "unpolished")])
    kept, withheld = apply_flank_carried_rule([just_inside, just_outside])
    assert [k.end for k in kept] == [14_000]
    assert [w.end for w in withheld] == [14_001]


def test_an_unmodelled_core_outside_the_flank_span_is_withheld():
    """The Debaryomyces case: the core fragment sits 8 kb beyond the block."""
    r = _result(FLANKS + [_ev("mfa1", "core_MAT", 19_000, 19_400, "unpolished")])
    kept, withheld = apply_flank_carried_rule([r])
    assert kept == []
    [w] = withheld
    assert w.withheld_reason == WITHHELD_FLANK_CARRIED


def test_one_core_hit_outside_is_enough_to_withhold():
    r = _result(FLANKS + [_ev("mfa1", "core_MAT", 5_000, 5_500, "unpolished"),
                          _ev("pra1", "core_MAT", 30_000, 30_500, "unpolished")])
    kept, withheld = apply_flank_carried_rule([r])
    assert kept == [] and len(withheld) == 1


def test_a_core_hit_on_another_contig_is_outside():
    r = _result(FLANKS + [_ev("mfa1", "core_MAT", 5_000, 5_500, "unpolished", contig="c2")])
    kept, withheld = apply_flank_carried_rule([r])
    assert kept == [] and len(withheld) == 1


def test_a_call_with_no_core_hit_at_all_is_withheld():
    """Nothing places a MAT gene between the flanks; there is no call to make."""
    kept, withheld = apply_flank_carried_rule([_result(FLANKS)])
    assert kept == [] and len(withheld) == 1


def test_a_call_with_no_modelled_flank_is_left_to_the_other_bars():
    """Not flank-carried: the modelled-gene bar and tiering already rule on it."""
    r = _result([_ev("flk1", "flanking_conserved", 1_000, 2_000, "unpolished"),
                 _ev("mfa1", "core_MAT", 5_000, 5_500, "unpolished")])
    kept, withheld = apply_flank_carried_rule([r])
    assert kept == [r] and withheld == []


ORDER = (
    "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
    "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
    "    genes:\n      - {name: mfa1, role: core_MAT}\n"
    "      - {name: flk1, role: flanking_conserved}\n"
    "      - {name: flk2, role: flanking_conserved}\n"
)


def _annotated(gene, role, start, end):
    return SearchHit(KEY, gene, role, "c1", start, end, "+", 91.5, "rec1",
                     "diamond_proteome", coverage=88.0)


def _raw(gene, start, end):
    return SearchHit(KEY, gene, "core_MAT", "c1", start, end, "+", 31.0, "rec1",
                     "tblastn_genome", coverage=12.0)


def _run(tmp_path, core_start, core_end):
    _write_order(tmp_path, ORDER)
    _write_record(tmp_path)
    hits = [_annotated("flk1", "flanking_conserved", 100, 200),
            _annotated("flk2", "flanking_conserved", 5_000, 5_100),
            _raw("mfa1", core_start, core_end)]
    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa",
        taxid=None, db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=lambda *a, **k: hits, search_localize=lambda *a, **k: [],
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )


def test_the_pipeline_keeps_an_inside_call_at_low(tmp_path):
    outcome = _run(tmp_path, 2_000, 2_300)
    [r] = outcome.results
    assert r.confidence == "low"
    assert r.locus_class == LOCUS_CLASS_PARTIAL
    assert r.idiomorph_unmodelled is True


def test_the_pipeline_withholds_an_outside_call_and_says_why(tmp_path):
    outcome = _run(tmp_path, 9_000, 9_300)
    assert outcome.results == []
    [w] = outcome.suppressed_loci
    assert w.withheld_reason == WITHHELD_FLANK_CARRIED
    assert outcome.suppressed_flank_carried == 1
    [reason] = [n.reason for n in outcome.not_detected if n.family_key == KEY]
    assert "flank" in reason

    path = tmp_path / "report.yaml"
    write_detection_report(outcome, path)
    doc = yaml.safe_load(path.read_text())
    assert doc["suppressed_flank_carried"] == 1
    assert doc["suppressed_loci"][0]["withheld_reason"] == WITHHELD_FLANK_CARRIED


def test_the_report_writes_the_idiomorph_unmodelled_flag(tmp_path):
    outcome = _run(tmp_path, 2_000, 2_300)
    path = tmp_path / "report.yaml"
    write_detection_report(outcome, path)
    doc = yaml.safe_load(path.read_text())
    assert doc["detected"][0]["idiomorph_unmodelled"] is True
