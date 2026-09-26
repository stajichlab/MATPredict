"""Withholding curated records, so recall can be measured rather than assumed.

A curated record is the QUERY the pipeline searches with, so running detection
on the genome it came from asks "did we find the thing we told it to look for".
`benchmark.py` has always reported `sensitivity=None` for every family because
nothing could withhold a record. These tests cover the piece that changes that.

RADIUS IS WHAT MAKES A RECALL NUMBER MEAN ANYTHING. Holding out one record
tests strain/assembly robustness and is nearly free to pass, since same-species
proteins run ~99% identical. Holding out an order tests whether an uncurated
clade can be called at all. The same "recall" figure means opposite things at
the two ends.
"""
from pathlib import Path

import pytest
import yaml

from MATPredict.detect.holdout import (
    Radius, _parse_lineage, load_record_taxa, records_to_withhold,
)
from MATPredict.detect.reference_fasta import build_reference_fasta

DB = Path("db")
pytestmark = pytest.mark.skipif(not DB.exists(), reason="needs the curated database")


def test_lineage_parses_into_ranks():
    r = _parse_lineage("k__Fungi;p__Ascomycota;o__Sordariales;g__Neurospora;s__Neurospora_crassa")
    assert r["o__"] == "Sordariales"
    assert r["s__"] == "Neurospora_crassa"


def test_a_malformed_lineage_yields_no_ranks_rather_than_guessing():
    # Silently widening a holdout because a lineage failed to parse would
    # OVERSTATE how hard the test was -- the direction that flatters the result.
    assert _parse_lineage("") == {}
    assert _parse_lineage("Fungi, Ascomycota") == {}


def test_the_subject_record_is_always_withheld():
    for radius in Radius:
        drop = records_to_withhold(DB, "5141_74-ors-a_MAT_MAT1-1", radius)
        assert "5141_74-ors-a_MAT_MAT1-1" in drop, radius


def test_radius_widens_monotonically():
    sizes = [len(records_to_withhold(DB, "746128_af293_MAT_MAT1-2", r))
             for r in (Radius.RECORD, Radius.SPECIES, Radius.GENUS,
                       Radius.FAMILY, Radius.ORDER)]
    assert sizes == sorted(sizes), sizes
    assert sizes[-1] > sizes[0], "an order holdout must be wider than a record one"


def test_species_radius_catches_the_other_idiomorph_of_the_same_species():
    # Neurospora crassa is curated twice, once per idiomorph. A species holdout
    # that kept the MAT1-2 record would leave a near-identical answer in place.
    drop = records_to_withhold(DB, "5141_74-ors-a_MAT_MAT1-1", Radius.SPECIES)
    assert {"5141_74-ors-a_MAT_MAT1-1", "5141_unknown-1_MAT_MAT1-2"} <= drop


def test_unknown_record_is_an_error_not_a_silent_no_op():
    with pytest.raises(KeyError):
        records_to_withhold(DB, "not_a_record", Radius.RECORD)


def test_withheld_records_leave_the_query_set(tmp_path):
    full = build_reference_fasta(DB, tmp_path / "full.faa")
    drop = records_to_withhold(DB, "5141_74-ors-a_MAT_MAT1-1", Radius.SPECIES)
    held = build_reference_fasta(DB, tmp_path / "held.faa", exclude_record_ids=drop)
    assert held.read_text().count(">") < full.read_text().count(">")
    for rid in drop:
        assert rid not in held.read_text()


def test_holding_out_everything_leaves_no_queries(tmp_path):
    all_ids = {t.record_id for t in load_record_taxa(DB)}
    out = build_reference_fasta(DB, tmp_path / "none.faa", exclude_record_ids=all_ids)
    assert out.read_text().count(">") == 0


def test_candidates_are_not_holdout_subjects():
    # db/candidates/ holds needs_review records that never reach the query set,
    # so they can neither be withheld nor leak an answer.
    for t in load_record_taxa(DB):
        assert t.path.relative_to(DB).parts[0] != "candidates"
