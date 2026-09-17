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


def test_write_detection_gff3_gene_parent_is_scoped_to_its_own_contig(tmp_path):
    """Finding B (part 1) regression: a fragmented locus's gene on the second
    contig must NOT carry a Parent pointing at a MAT_locus feature declared
    only on the first contig -- that is not valid/clean GFF3 for a
    multi-contig feature set. Each segment gets its own MAT_locus feature,
    scoped to its own contig, and each gene's Parent points at the segment
    feature sharing its contig."""
    fragmented = DetectionResult(
        family_key=KEY, contig="c1", start=100, end=400,
        confidence="medium", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["mfa1", "pra1"], genes_missing=[], fragmented=True,
        segments=[LocusSegment("c1", 100, 200, 99), LocusSegment("c2", 300, 400, 299)],
        gene_evidence=[
            GeneEvidence("mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, None,
                         "5270_521_aLocus_a1", "diamond_proteome"),
            GeneEvidence("pra1", "core_MAT", "c2", 300, 400, "+", 95.0, None,
                         "5270_521_aLocus_a1", "diamond_proteome"),
        ],
    )
    out = tmp_path / "frag.gff3"
    write_detection_gff3(DetectionOutcome(results=[fragmented]), out)
    lines = out.read_text().splitlines()

    feature_lines = [line for line in lines if not line.startswith("#")]
    locus_lines = [line for line in feature_lines if "\tMAT_locus\t" in line]
    gene_lines = [line for line in feature_lines if "\tgene\t" in line]

    def _attr(line: str, key: str) -> str:
        attrs = dict(a.split("=", 1) for a in line.split("\t")[8].split(";"))
        return attrs[key]

    # one MAT_locus feature per segment, each declared on its own contig
    assert len(locus_lines) == 2
    locus_by_contig = {line.split("\t")[0]: line for line in locus_lines}
    c1_locus_id = _attr(locus_by_contig["c1"], "ID")
    c2_locus_id = _attr(locus_by_contig["c2"], "ID")
    assert c1_locus_id != c2_locus_id

    # every gene's Parent is on ITS OWN contig, never the other segment's contig
    mfa1_line = next(line for line in gene_lines if "Name=mfa1" in line)
    pra1_line = next(line for line in gene_lines if "Name=pra1" in line)
    assert mfa1_line.split("\t")[0] == "c1"
    assert _attr(mfa1_line, "Parent") == c1_locus_id
    assert pra1_line.split("\t")[0] == "c2"
    assert _attr(pra1_line, "Parent") == c2_locus_id


def test_write_detection_gff3_deduplicates_sequence_region_across_results(tmp_path):
    """Finding B (part 2) regression: two separate DetectionResults that both
    reference contig c1 must not each emit their own ##sequence-region c1
    pragma -- GFF3 tooling expects at most one per seqid."""
    key_b = FamilyKey("Basidiomycota", "bLocus")
    result_a = DetectionResult(
        family_key=KEY, contig="c1", start=100, end=200,
        confidence="high", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["mfa1"], genes_missing=[], fragmented=False,
    )
    result_b = DetectionResult(
        family_key=key_b, contig="c1", start=500, end=600,
        confidence="high", idiomorph="undetermined", ambiguous_with=[],
        genes_found=["bE"], genes_missing=[], fragmented=False,
    )
    out = tmp_path / "dup.gff3"
    write_detection_gff3(DetectionOutcome(results=[result_a, result_b]), out)
    text = out.read_text()
    assert text.count("##sequence-region c1") == 1
    # the deduplicated pragma widens to cover both results' extents
    assert "##sequence-region c1 100 600" in text


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
