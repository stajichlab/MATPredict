"""The losing half of a resolved pair is chosen from the events that agree with the locus call.

Xylaria bambusicola (GCA_033171785.1): MAT1-2-1 beat MAT1-1-1, MAT1-1-2 and MAT1-1-3 where its model sits, but a weak
hit 15 kb away was won by MAT1-1-5. Keyed by gene name, that far event marked the polished MAT1-2-1 model "lost to
MAT1-1-5", and the GFF and figure (which omit losers) dropped the locus's strongest MAT gene."""
from __future__ import annotations

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.idiomorph import IdiomorphResolution, lost_to_for_call

FAMILY = Family(FamilyKey("Ascomycota", "MAT"), "enum", ["MAT1-1", "MAT1-2"], None, [
    {"name": "MAT1-1-1", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-1"]},
    {"name": "MAT1-1-3", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-1"]},
    {"name": "MAT1-1-5", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-1"]},
    {"name": "MAT1-2-1", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-2"]},
    {"name": "APN2", "role": "flanking_conserved"},
], [4890])


def _event(winner, loser):
    return IdiomorphResolution(contig="c", winner=winner, loser=loser, winner_identity=40.0, loser_identity=30.0, overlap_fraction=1.0)


XYLARIA = [_event("MAT1-1-1", "MAT1-2-1"), _event("MAT1-2-1", "MAT1-1-1"), _event("MAT1-2-1", "MAT1-1-3"),
           _event("MAT1-1-3", "MAT1-2-1"), _event("MAT1-1-5", "MAT1-2-1")]


def test_the_called_idiomorphs_gene_is_never_marked_lost():
    lost = lost_to_for_call(FAMILY, "MAT1-2", XYLARIA)
    assert "MAT1-2-1" not in lost
    assert lost == {"MAT1-1-1": "MAT1-2-1", "MAT1-1-3": "MAT1-2-1"}


def test_the_other_calls_events_are_the_ones_kept():
    lost = lost_to_for_call(FAMILY, "MAT1-1", XYLARIA)
    assert lost == {"MAT1-2-1": "MAT1-1-5"}          # the last agreeing event, as before


def test_an_undetermined_call_keeps_every_event_as_before():
    legacy = {"MAT1-2-1": "MAT1-1-5", "MAT1-1-1": "MAT1-2-1", "MAT1-1-3": "MAT1-2-1"}
    assert lost_to_for_call(FAMILY, "undetermined", XYLARIA) == legacy
    assert lost_to_for_call(FAMILY, None, XYLARIA) == legacy


def test_a_gene_with_no_idiomorph_does_not_block_an_event():
    # APN2 is a flank gene: it excludes nothing, so an event naming it is kept whatever the call.
    assert lost_to_for_call(FAMILY, "MAT1-2", [_event("MAT1-2-1", "APN2")]) == {"APN2": "MAT1-2-1"}
