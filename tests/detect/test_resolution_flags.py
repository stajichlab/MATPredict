"""`resolution_flags`: a called locus says when it may be a fragment; calls are unchanged (#256)."""
from __future__ import annotations

import yaml

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import DetectionOutcome, DetectionResult, GeneEvidence, LocusSegment
from MATPredict.detect.report import write_detection_report
from MATPredict.detect.resolution_flags import resolution_flags

KEY = FamilyKey("Ascomycota", "MATog")


def _gene(name, start, end, coverage=90.0, role="core_MAT", contig="c1", exons=None, lost_to=None):
    return GeneEvidence(name, role, contig, start, end, "+", 60.0, coverage, "rec", "exonerate_refine", exons=exons, lost_to=lost_to)


def test_an_intact_mid_contig_core_gene_has_no_flag():
    assert resolution_flags([_gene("MATa1", 5000, 5600)], {"c1": 100000}) == {}


def test_a_core_gene_at_a_contig_end_is_flagged_either_end():
    assert resolution_flags([_gene("MATa1", 5000, 5600)], {"c1": 5900}) == {"contig_edge": ["MATa1"]}
    assert resolution_flags([_gene("MATa1", 100, 700)], {"c1": 100000}) == {"contig_edge": ["MATa1"]}


def test_unknown_contig_length_gives_no_edge_flag():
    assert resolution_flags([_gene("MATa1", 100, 700)], None) == {}


def test_low_coverage_and_short_models_are_flagged_and_flanks_are_not():
    g = [_gene("MATa1", 5000, 5600, coverage=30.0), _gene("MATalpha2", 8000, 8150), _gene("sla2", 100, 160, coverage=5.0, role="flanking_conserved")]
    assert resolution_flags(g, {"c1": 100000}) == {"partial_core": ["MATa1"], "short_core": ["MATalpha2"]}


def test_exons_set_the_model_length_and_the_losing_half_of_a_pair_is_skipped():
    spliced = _gene("MATalpha2", 8000, 9000, exons=[(8000, 8100), (8900, 9000)])        # 202 bp = 67 aa
    assert resolution_flags([spliced], {"c1": 100000}) == {"short_core": ["MATalpha2"]}
    assert resolution_flags([_gene("MATa1", 100, 150, lost_to="MATalpha1")], {"c1": 100000}) == {}


def test_the_report_carries_the_flags_and_the_call_is_untouched(tmp_path):
    genes = [_gene("MATa1", 1, 170, coverage=30.0)]
    r = DetectionResult(family_key=KEY, contig="c1", start=1, end=170, confidence="high", idiomorph="MATa", ambiguous_with=[], genes_found=["MATa1"],
                        genes_missing=[], fragmented=False, segments=[LocusSegment("c1", 1, 170, contig_edge_distance=0)], gene_evidence=genes, reference_records=["rec"])
    out = tmp_path / "r.yaml"
    write_detection_report(DetectionOutcome(results=[r], not_detected=[], families_attempted=[KEY], contig_lengths={"c1": 1000}), out)
    d = yaml.safe_load(out.read_text())["detected"][0]
    assert d["resolution_flags"] == ["contig_edge", "partial_core", "short_core"]
    assert d["resolution_flag_genes"]["contig_edge"] == ["MATa1"] and d["confidence"] == "high" and d["idiomorph"] == "MATa"
