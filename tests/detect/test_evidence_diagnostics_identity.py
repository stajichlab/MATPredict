# tests/detect/test_evidence_diagnostics_identity.py
"""Diagnostics rows say which run and which genome produced them.

This file is the calibration corpus -- the thing the provisional 0.5 overlap
threshold and any future `min_identity` are supposed to be derived from. It
could not serve that purpose: a row carried family, contig, cluster span,
counts, roles and identity, but nothing identifying the GENOME. Concatenate a
batch's rows and there is no way to tell one organism's clusters from
another's, so no per-genome statistic can be computed from it at all.

The writer also appends and never truncates, which is correct for a batch
sharing one file but means a re-run silently doubles the corpus. That is the
real explanation for the "60 rows = 30 unique pairs" noted during validation:
across all 67 files of the real sweep there are 2335 rows and ZERO duplicates,
so nothing was ever double-emitted. A run id makes an accidental second run
visible and de-duplicable instead of invisible.
"""
from __future__ import annotations

import json

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.idiomorph import IdiomorphResolution
from MATPredict.detect.pipeline import (
    _write_evidence_diagnostics,
    _write_idiomorph_diagnostics,
)
from MATPredict.detect.search import SearchHit

FAM = Family(
    FamilyKey("Mucoromycota", "MAT"), "enum", ["Plus", "Minus"], None,
    [
        {"name": "sexP", "role": "core_MAT", "present_in_idiomorphs": ["Plus"]},
        {"name": "sexM", "role": "core_MAT", "present_in_idiomorphs": ["Minus"]},
    ],
    [4827],
)


def _cluster():
    hit = SearchHit(
        FAM.key, "sexP", "core_MAT", "c1", 1, 100, "+", 47.3, "rec1",
        "diamond_proteome", coverage=23.6,
    )
    return GeneCluster("c1", 1, 100, [hit])


def _rows(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def test_an_evidence_row_names_its_run_and_its_genome(tmp_path):
    out = tmp_path / "diag.jsonl"
    _write_evidence_diagnostics(
        out, _cluster(), FAM, admitted=True, run_id="r1", genome_id="Absidia_x",
    )
    row = _rows(out)[0]
    assert row["run_id"] == "r1"
    assert row["genome_id"] == "Absidia_x"
    assert row["kind"] == "evidence"


def test_two_runs_appending_to_one_file_stay_distinguishable(tmp_path):
    # Append is deliberate -- a batch writes many genomes to one file -- so
    # the run id is what makes an accidental re-run detectable rather than
    # silently doubling the corpus.
    out = tmp_path / "diag.jsonl"
    _write_evidence_diagnostics(out, _cluster(), FAM, True, "r1", "Absidia_x")
    _write_evidence_diagnostics(out, _cluster(), FAM, True, "r2", "Absidia_x")
    assert [r["run_id"] for r in _rows(out)] == ["r1", "r2"]


def test_a_resolution_event_records_everything_a_recalibration_needs(tmp_path):
    out = tmp_path / "diag.jsonl"
    _write_idiomorph_diagnostics(
        out,
        IdiomorphResolution(
            contig="c1", winner="sexM", loser="sexP",
            winner_identity=30.60, loser_identity=28.26,
            overlap_fraction=0.97, winner_coverage=None, loser_coverage=23.6,
        ),
        FAM, run_id="r1", genome_id="Cunninghamella_x",
    )
    row = _rows(out)[0]
    assert row["kind"] == "idiomorph_resolution"
    assert row["run_id"] == "r1"
    assert row["genome_id"] == "Cunninghamella_x"
    assert row["family"] == "Mucoromycota:MAT"
    assert row["winner"] == "sexM"
    assert row["loser"] == "sexP"
    assert row["winner_identity"] == 30.60
    assert row["loser_identity"] == 28.26
    assert abs(row["margin"] - 2.34) < 1e-9
    assert row["overlap_fraction"] == 0.97
    # Coverage on BOTH members, because whether coverage discriminates better
    # than identity is the open question this corpus exists to answer.
    assert row["winner_coverage"] is None
    assert row["loser_coverage"] == 23.6


def test_a_diagnostics_write_failure_never_aborts_a_run(tmp_path):
    # Diagnostics are best-effort. Losing the corpus is bad; losing a
    # detection run because a log path was wrong is worse.
    missing = tmp_path / "no_such_dir" / "diag.jsonl"
    _write_evidence_diagnostics(missing, _cluster(), FAM, True, "r1", "g1")
    _write_idiomorph_diagnostics(
        missing,
        IdiomorphResolution("c1", "sexM", "sexP", 30.6, 28.26, 1.0),
        FAM, run_id="r1", genome_id="g1",
    )
