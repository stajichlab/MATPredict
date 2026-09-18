from __future__ import annotations

from pathlib import Path

import yaml

from MATPredict.detect.rollout_aggregate import (
    Anomaly,
    GenomeReportError,
    NotDetectedEntry,
    aggregate_reports,
    write_rollout_summary,
)

# Real-shaped fixtures matching `write_detection_report`'s exact output shape
# (`_result_doc`/`write_detection_report` in report.py): `families_attempted`,
# `detected: [{family, confidence, ...}]`, `not_detected: [{family, reason, ...}]`.

HIGH_CONFIDENCE_DOC = {
    "families_attempted": ["Basidiomycota:aLocus"],
    "detected": [
        {
            "family": "Basidiomycota:aLocus",
            "contig": "c1",
            "start": 100,
            "end": 6443,
            "confidence": "high",
            "idiomorph": "undetermined",
            "ambiguous_with": [],
            "genes_found": ["pra1"],
            "genes_missing": [],
            "genes_not_searchable": [],
            "fragmented": False,
            "reference_records": ["5270_521_aLocus_a1"],
            "segments": [{"contig": "c1", "start": 100, "end": 6443, "contig_edge_distance": 99}],
            "gene_evidence": [
                {
                    "gene": "pra1", "role": "core_MAT", "contig": "c1", "start": 1000, "end": 2000,
                    "strand": "+", "identity": 92.5, "coverage": 87.0,
                    "reference_record": "5270_521_aLocus_a1", "method": "diamond_proteome",
                    "status": "found", "alternate_model": None,
                }
            ],
        }
    ],
    "not_detected": [],
}

NOT_DETECTED_DOC = {
    "families_attempted": ["Basidiomycota:aLocus", "Basidiomycota:bLocus"],
    "detected": [],
    "not_detected": [
        {
            "family": "Basidiomycota:bLocus",
            "reason": "best cluster matched 0.25 of this family's expected genes, "
                      "below the ambiguity floor of 0.50",
            "best_fraction_found": 0.25,
            "genes_found": ["bE"],
            "genes_missing": ["bW"],
            "genes_not_searchable": [],
        }
    ],
}

LOW_TIER_DOC = {
    "families_attempted": ["Basidiomycota:aLocus"],
    "detected": [
        {
            "family": "Basidiomycota:aLocus",
            "contig": "c2",
            "start": 200,
            "end": 5000,
            "confidence": "low",
            "idiomorph": "undetermined",
            "ambiguous_with": [],
            "genes_found": ["pra1"],
            "genes_missing": ["rba1"],
            "genes_not_searchable": [],
            "fragmented": False,
            "reference_records": ["5270_521_aLocus_a1"],
            "segments": [{"contig": "c2", "start": 200, "end": 5000, "contig_edge_distance": 10}],
            "gene_evidence": [],
        }
    ],
    "not_detected": [],
}


def _write(tmp_path: Path, genome_id: str, doc: dict) -> Path:
    genome_dir = tmp_path / genome_id
    genome_dir.mkdir()
    report_path = genome_dir / "detection_report.yaml"
    report_path.write_text(yaml.safe_dump(doc, sort_keys=False))
    return report_path


def test_aggregate_reports_tallies_confidence_and_surfaces_not_detected_reasons(tmp_path):
    p1 = _write(tmp_path, "1_GCA_000000001.1", HIGH_CONFIDENCE_DOC)
    p2 = _write(tmp_path, "2_GCA_000000002.1", NOT_DETECTED_DOC)
    p3 = _write(tmp_path, "3_GCA_000000003.1", LOW_TIER_DOC)

    summary = aggregate_reports([p1, p2, p3])

    assert summary.total_genomes == 3
    assert summary.confidence_tally == {
        "Basidiomycota:aLocus": {"high": 1, "low": 1},
    }
    assert summary.not_detected == [
        NotDetectedEntry(
            genome="2_GCA_000000002.1",
            family="Basidiomycota:bLocus",
            reason="best cluster matched 0.25 of this family's expected genes, "
                   "below the ambiguity floor of 0.50",
        )
    ]
    assert summary.genome_errors == []


def test_aggregate_reports_treats_missing_report_as_no_result_not_a_crash(tmp_path):
    p1 = _write(tmp_path, "1_GCA_000000001.1", HIGH_CONFIDENCE_DOC)
    # Task 3's known gap: run_batch creates the genome directory before its
    # pipeline call, so a failed genome leaves an empty directory with no
    # detection_report.yaml inside it.
    empty_genome_dir = tmp_path / "4_GCA_000000004.1"
    empty_genome_dir.mkdir()
    missing_report = empty_genome_dir / "detection_report.yaml"

    summary = aggregate_reports([p1, missing_report])

    assert summary.total_genomes == 2
    assert summary.confidence_tally == {"Basidiomycota:aLocus": {"high": 1}}
    assert len(summary.genome_errors) == 1
    assert summary.genome_errors[0].genome == "4_GCA_000000004.1"


def test_aggregate_reports_flags_anomaly_when_relative_misses_a_detected_family(tmp_path):
    detected_doc = {
        "families_attempted": ["Basidiomycota:aLocus"],
        "detected": [dict(HIGH_CONFIDENCE_DOC["detected"][0])],
        "not_detected": [],
    }
    missed_doc = {
        "families_attempted": ["Basidiomycota:aLocus"],
        "detected": [],
        "not_detected": [
            {
                "family": "Basidiomycota:aLocus",
                "reason": "no candidate cluster found",
                "best_fraction_found": 0.0,
                "genes_found": [],
                "genes_missing": ["pra1"],
                "genes_not_searchable": [],
            }
        ],
    }
    p1 = _write(tmp_path, "111_GCA_AAA.1", detected_doc)
    p2 = _write(tmp_path, "222_GCA_BBB.1", missed_doc)

    def fake_lineage(taxid: int) -> str:
        # Both taxids resolve to the same order, so the second genome's miss
        # is anomalous relative to the first genome's detection.
        return "k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__X;g__Y;s__Z"

    summary = aggregate_reports([p1, p2], lineage_resolver=fake_lineage)

    assert summary.anomalies == [
        Anomaly(
            genome="222_GCA_BBB.1",
            family="Basidiomycota:aLocus",
            taxonomic_group="order:Agaricales",
            detected_in=["111_GCA_AAA.1"],
        )
    ]


def test_aggregate_reports_no_anomaly_when_lineage_resolver_fails(tmp_path):
    detected_doc = {
        "families_attempted": ["Basidiomycota:aLocus"],
        "detected": [dict(HIGH_CONFIDENCE_DOC["detected"][0])],
        "not_detected": [],
    }
    missed_doc = {
        "families_attempted": ["Basidiomycota:aLocus"],
        "detected": [],
        "not_detected": [
            {
                "family": "Basidiomycota:aLocus",
                "reason": "no candidate cluster found",
                "best_fraction_found": 0.0,
                "genes_found": [],
                "genes_missing": ["pra1"],
                "genes_not_searchable": [],
            }
        ],
    }
    p1 = _write(tmp_path, "111_GCA_AAA.1", detected_doc)
    p2 = _write(tmp_path, "222_GCA_BBB.1", missed_doc)

    def unresolvable(taxid: int):
        raise RuntimeError("taxonkit not installed")

    summary = aggregate_reports([p1, p2], lineage_resolver=unresolvable)

    assert summary.anomalies == []


def test_write_rollout_summary_round_trips_as_yaml(tmp_path):
    p1 = _write(tmp_path, "1_GCA_000000001.1", HIGH_CONFIDENCE_DOC)
    summary = aggregate_reports([p1])

    out_path = tmp_path / "rollout_summary.yaml"
    write_rollout_summary(summary, out_path)

    doc = yaml.safe_load(out_path.read_text())
    assert doc["total_genomes"] == 1
    assert doc["confidence_tally"] == {"Basidiomycota:aLocus": {"high": 1}}
