from __future__ import annotations
import yaml

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import DetectionResult
from MATPredict.detect.report import write_detection_gff3, write_detection_report

RESULT = DetectionResult(
    family_key=FamilyKey("Basidiomycota", "aLocus"), contig="c1", start=100, end=6443,
    confidence="high", idiomorph="undetermined", ambiguous_with=[],
    genes_found=["mfa1", "pra1"], genes_missing=[], fragmented=False,
)


def test_write_detection_gff3(tmp_path):
    out = tmp_path / "out.gff3"
    write_detection_gff3([RESULT], out)
    text = out.read_text()
    assert text.startswith("##gff-version 3")
    assert "c1\tMATPredict\tMAT_locus\t100\t6443" in text


def test_write_detection_report(tmp_path):
    out = tmp_path / "report.yaml"
    write_detection_report([RESULT], out)
    doc = yaml.safe_load(out.read_text())
    assert doc[0]["family"] == "Basidiomycota:aLocus"
    assert doc[0]["confidence"] == "high"
    assert doc[0]["genes_missing"] == []
