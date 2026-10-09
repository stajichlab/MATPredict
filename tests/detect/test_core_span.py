"""`core_span`: where a locus's own genes lie, reported next to the cluster span (report only; calls are unchanged).

Why: the cluster span is built from every hit of every family, so one weak hit from another family's short query can
stretch a locus by 20 kb (T48-F, analysis/2026-10-08_cinerea-b43-trace.md). core_span shows how far the span goes beyond the
locus's own genes. Seam: the detect report writers, fed a DetectionOutcome.
"""
from __future__ import annotations

import yaml

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import DetectionOutcome, DetectionResult, GeneEvidence, LocusSegment
from MATPredict.detect.report import write_detection_gff3, write_detection_report

KEY = FamilyKey("Basidiomycota", "HD")


def _gene(name, contig, start, end):
    return GeneEvidence(name, "core_MAT", contig, start, end, "+", 80.0, 90.0, "rec_HD_A1", "exonerate_refine")


def _result(genes, start=100, end=6443):
    return DetectionResult(
        family_key=KEY, contig="c1", start=start, end=end, confidence="high", idiomorph="undetermined",
        ambiguous_with=[], genes_found=[g.gene_name for g in genes], genes_missing=[], fragmented=False,
        segments=[LocusSegment("c1", start, end, contig_edge_distance=99)], gene_evidence=genes,
        reference_records=["rec_HD_A1"],
    )


def _outcome(result):
    return DetectionOutcome(results=[result], not_detected=[], families_attempted=[KEY])


def test_report_gives_the_extent_of_the_locus_own_genes_on_its_contig(tmp_path):
    # cluster span 100-6443 (6,344 bp); own genes at 1000-2000 and 3000-4000; one gene on another contig is ignored
    genes = [_gene("HD1", "c1", 1000, 2000), _gene("HD2", "c1", 3000, 4000), _gene("MIP1", "c2", 50, 90)]
    out = tmp_path / "r.yaml"
    write_detection_report(_outcome(_result(genes)), out)
    locus = yaml.safe_load(out.read_text())["detected"][0]
    assert locus["core_span"] == {"start": 1000, "end": 4000, "beyond_core_bp": 3343}


def test_a_span_that_equals_its_genes_has_nothing_beyond_the_core(tmp_path):
    genes = [_gene("HD1", "c1", 100, 3000), _gene("HD2", "c1", 3500, 6443)]
    out = tmp_path / "r.yaml"
    write_detection_report(_outcome(_result(genes)), out)
    assert yaml.safe_load(out.read_text())["detected"][0]["core_span"] == {"start": 100, "end": 6443, "beyond_core_bp": 0}


def test_a_locus_without_gene_models_has_no_core_span(tmp_path):
    out = tmp_path / "r.yaml"
    write_detection_report(_outcome(_result([])), out)
    assert yaml.safe_load(out.read_text())["detected"][0]["core_span"] is None


def test_gff3_locus_line_carries_the_core_extent(tmp_path):
    genes = [_gene("HD1", "c1", 1000, 2000), _gene("HD2", "c1", 3000, 4000)]
    out = tmp_path / "o.gff3"
    write_detection_gff3(_outcome(_result(genes)), out)
    locus = next(l for l in out.read_text().splitlines() if "\tMAT_locus\t" in l)
    attrs = dict(kv.split("=", 1) for kv in locus.split("\t")[8].split(";") if "=" in kv)
    assert (locus.split("\t")[3], locus.split("\t")[4]) == ("100", "6443")          # the cluster span is unchanged
    assert (attrs["core_start"], attrs["core_end"], attrs["beyond_core_bp"]) == ("1000", "4000", "3343")


def test_gff3_locus_line_has_no_core_attributes_without_gene_models(tmp_path):
    out = tmp_path / "o.gff3"
    write_detection_gff3(_outcome(_result([])), out)
    locus = next(l for l in out.read_text().splitlines() if "\tMAT_locus\t" in l)
    assert "core_start" not in locus and "beyond_core_bp" not in locus
