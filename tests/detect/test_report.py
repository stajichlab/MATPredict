from __future__ import annotations
import yaml

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import (
    DetectionOutcome,
    DetectionResult,
    GeneEvidence,
    LocusSegment,
    NotDetectedFamily,
)
from MATPredict.detect.report import write_detection_gff3, write_detection_report

KEY = FamilyKey("Basidiomycota", "aLocus")

RESULT = DetectionResult(
    family_key=KEY, contig="c1", start=100, end=6443,
    confidence="high", idiomorph="undetermined", ambiguous_with=[],
    genes_found=["pra1"], genes_missing=["rba1"], fragmented=False,
    genes_not_searchable=["mfa1"],
    segments=[LocusSegment("c1", 100, 6443, contig_edge_distance=99)],
    gene_evidence=[
        GeneEvidence("pra1", "core_MAT", "c1", 1000, 2000, "+", 92.5, 87.0,
                     "5270_521_aLocus_a1", "diamond_proteome"),
    ],
    reference_records=["5270_521_aLocus_a1"],
)

OUTCOME = DetectionOutcome(
    results=[RESULT],
    not_detected=[
        NotDetectedFamily(
            family_key=FamilyKey("Basidiomycota", "bLocus"),
            reason="best cluster matched 0.25 of this family's expected genes, "
                   "below the ambiguity floor of 0.50",
            best_fraction_found=0.25,
            genes_found=["bE"],
            genes_missing=["bW"],
        )
    ],
    families_attempted=[KEY, FamilyKey("Basidiomycota", "bLocus")],
)


def test_write_detection_gff3_emits_locus_and_gene_features(tmp_path):
    out = tmp_path / "out.gff3"
    write_detection_gff3(OUTCOME, out)
    text = out.read_text()
    assert text.startswith("##gff-version 3")
    assert "##sequence-region c1 100 6443" in text
    assert "c1\tMATPredict\tMAT_locus\t100\t6443" in text
    # per-gene FEATURE lines, not just one locus-region line
    gene_line = next(line for line in text.splitlines() if "\tgene\t" in line and "Name=pra1" in line)
    assert "\t1000\t2000\t" in gene_line
    assert "role=core_MAT" in gene_line
    assert "present=true" in gene_line
    assert "identity=92.5" in gene_line
    assert "coverage=87.0" in gene_line
    assert "reference_record=5270_521_aLocus_a1" in gene_line
    # absent and not-searchable genes are explicit, with distinct semantics
    assert "Name=rba1;present=false" in text
    assert "Name=mfa1;present=false;not_searchable=true" in text


def test_write_detection_gff3_emits_one_sequence_region_per_segment(tmp_path):
    fragmented = DetectionResult(
        family_key=KEY, contig="c1", start=100, end=400,
        confidence="medium", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["mfa1", "pra1"], genes_missing=[], fragmented=True,
        segments=[LocusSegment("c1", 100, 200, 99), LocusSegment("c2", 300, 400, 299)],
    )
    out = tmp_path / "frag.gff3"
    write_detection_gff3(DetectionOutcome(results=[fragmented]), out)
    text = out.read_text()
    assert "##sequence-region c1 100 200" in text
    assert "##sequence-region c2 300 400" in text
    assert "fragmented=true" in text


def test_write_detection_report(tmp_path):
    out = tmp_path / "report.yaml"
    write_detection_report(OUTCOME, out)
    doc = yaml.safe_load(out.read_text())
    detected = doc["detected"][0]
    assert detected["family"] == "Basidiomycota:aLocus"
    assert detected["confidence"] == "high"
    assert detected["genes_missing"] == ["rba1"]
    assert detected["genes_not_searchable"] == ["mfa1"]
    assert detected["reference_records"] == ["5270_521_aLocus_a1"]
    assert detected["segments"][0]["contig_edge_distance"] == 99
    evidence = detected["gene_evidence"][0]
    assert evidence["gene"] == "pra1"
    assert evidence["identity"] == 92.5
    assert evidence["coverage"] == 87.0
    assert evidence["reference_record"] == "5270_521_aLocus_a1"


def test_write_detection_report_lists_not_detected_families(tmp_path):
    """Sub-floor families must appear with a reason, never be silently dropped."""
    out = tmp_path / "report.yaml"
    write_detection_report(OUTCOME, out)
    doc = yaml.safe_load(out.read_text())
    assert doc["families_attempted"] == ["Basidiomycota:aLocus", "Basidiomycota:bLocus"]
    assert doc["not_detected"][0]["family"] == "Basidiomycota:bLocus"
    assert "below the ambiguity floor" in doc["not_detected"][0]["reason"]
    assert doc["not_detected"][0]["best_fraction_found"] == 0.25
