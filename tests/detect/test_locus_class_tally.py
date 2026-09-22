"""`locus_class_tally` -- the class breakdown beside `confidence_tally`.

Curator's ruling, 2026-09-21. `partial_locus` calls COUNT as detections in
recall, so the aggregate must show what fraction of a rollout's output is
partial rather than merging it into one locus count. Without this, a sweep
over 283 or 3,174 genomes reports a single number that could be mostly
partial, which is precisely the inflated statistic this project does not want.

The field is already written into every per-genome report, so this is a tally,
not new data collection.
"""
from pathlib import Path

import yaml

from MATPredict.detect.rollout_aggregate import aggregate_reports


def _report(tmp_path: Path, genome: str, detected: list[dict]) -> Path:
    path = tmp_path / f"{genome}.yaml"
    path.write_text(yaml.safe_dump({"genome": genome, "detected": detected}))
    return path


def test_locus_classes_are_tallied_per_family(tmp_path):
    paths = [
        _report(tmp_path, "g1", [
            {"family": "Mucoromycota/MAT", "confidence": "high", "locus_class": "mat_locus"},
            {"family": "Mucoromycota/MAT", "confidence": "medium", "locus_class": "partial_locus"},
        ]),
        _report(tmp_path, "g2", [
            {"family": "Mucoromycota/MAT", "confidence": "medium", "locus_class": "partial_locus"},
            {"family": "Basidiomycota/HD", "confidence": "high", "locus_class": "mat_locus"},
        ]),
    ]
    summary = aggregate_reports(paths)
    assert summary.locus_class_tally == {
        "Mucoromycota/MAT": {"mat_locus": 1, "partial_locus": 2},
        "Basidiomycota/HD": {"mat_locus": 1},
    }


def test_the_tally_reaches_the_serialized_document(tmp_path):
    paths = [_report(tmp_path, "g1", [
        {"family": "Mucoromycota/MAT", "confidence": "high", "locus_class": "mat_locus"},
    ])]
    doc = aggregate_reports(paths).to_doc()
    assert doc["locus_class_tally"] == {"Mucoromycota/MAT": {"mat_locus": 1}}


def test_a_result_with_no_locus_class_is_still_counted_in_confidence(tmp_path):
    """Older reports predate the field. They must not vanish from the
    confidence tally just because the class is missing."""
    paths = [_report(tmp_path, "g1", [
        {"family": "Mucoromycota/MAT", "confidence": "high"},
    ])]
    summary = aggregate_reports(paths)
    assert summary.confidence_tally == {"Mucoromycota/MAT": {"high": 1}}
    assert summary.locus_class_tally == {}


def test_partial_and_confirmed_are_never_merged(tmp_path):
    """The whole point: a reader must be able to see the split."""
    detected = [
        {"family": "F", "confidence": "medium", "locus_class": "partial_locus"}
        for _ in range(9)
    ] + [{"family": "F", "confidence": "high", "locus_class": "mat_locus"}]
    summary = aggregate_reports([_report(tmp_path, "g1", detected)])
    assert summary.locus_class_tally["F"] == {"partial_locus": 9, "mat_locus": 1}
    assert sum(summary.locus_class_tally["F"].values()) == 10
