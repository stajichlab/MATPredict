"""Withheld loci carry their classifier verdict (review finding F5, part).

Without it a withheld locus's idiomorph margin cannot be audited from the
report (results/2026-09-28_fable_review/).
"""
from dataclasses import replace

import yaml

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import MODELLED_GENE_BAR, DetectionOutcome, DetectionResult
from MATPredict.detect.report import write_detection_report


def _withheld(block):
    r = DetectionResult(FamilyKey("Mucoromycota", "MAT"), "c1", 1, 100, "medium",
                        "undetermined", [], ["sexM"], [], False)
    return replace(r, withheld_reason=MODELLED_GENE_BAR, idiomorph_classifier=block)


def test_a_withheld_locus_carries_its_classifier_block(tmp_path):
    block = {"verdict": "undetermined", "margin": 19.5, "input": "hsp_fragment"}
    out = tmp_path / "r.yaml"
    write_detection_report(DetectionOutcome(results=[], suppressed_loci=[_withheld(block)]), out)
    (row,) = yaml.safe_load(out.read_text())["suppressed_loci"]
    assert row["idiomorph_classifier"] == block


def test_a_withheld_locus_without_a_classifier_writes_null(tmp_path):
    out = tmp_path / "r.yaml"
    write_detection_report(DetectionOutcome(results=[], suppressed_loci=[_withheld(None)]), out)
    (row,) = yaml.safe_load(out.read_text())["suppressed_loci"]
    assert row["idiomorph_classifier"] is None
