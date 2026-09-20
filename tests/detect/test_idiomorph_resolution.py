# tests/detect/test_idiomorph_resolution.py
"""Resolving two mutually exclusive idiomorph genes that hit the SAME locus gene.

`sexM` and `sexP` share an HMG box, so a single real locus gene is hit by both
curated references. Measured on all 23 ground-truth Mucoromycota genomes this
happens every time, at 92-100% coordinate overlap, and it makes the idiomorph
uncallable: two idiomorphs are always indicated, so `assign_idiomorph` returns
`undetermined` in 23/23 and `expected_genes_for_idiomorph` never narrows.

The discriminator is coordinate overlap, not the mere co-occurrence of the two
gene names. Both idiomorphs genuinely present at one locus is real biology that
this project's curated records already document; what is NOT real is one
predicted protein being counted as two genes. Two hits to the same protein are
one gene; two hits to distinct, non-overlapping proteins are two.

The loser is marked superseded rather than dropped, so it stays in the report
as evidence of the ambiguity while every downstream consumer -- the evidence
floor, scoring, idiomorph assignment -- sees one gene.
"""
from __future__ import annotations

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.idiomorph import resolve_idiomorph_overlaps
from MATPredict.detect.search import SearchHit

FAM = Family(
    FamilyKey("Mucoromycota", "MAT"), "enum", ["Plus", "Minus"], None,
    [
        {"name": "tptA", "role": "flanking_conserved"},
        {"name": "sexP", "role": "core_MAT", "present_in_idiomorphs": ["Plus"]},
        {"name": "sexM", "role": "core_MAT", "present_in_idiomorphs": ["Minus"]},
        {"name": "rnhA", "role": "flanking_conserved"},
    ],
    [4827],
)


def _hit(gene, start, end, identity, role="core_MAT", coverage=None, contig="c1"):
    return SearchHit(
        FAM.key, gene, role, contig, start, end, "+", identity,
        f"rec_{gene}", "diamond_proteome", coverage=coverage,
    )


def _by_gene(hits):
    return {h.gene_name: h for h in hits}


def test_the_lower_identity_member_of_an_overlapping_pair_is_superseded():
    # The real Absidia cuneospora numbers: sexP 16939-17547 at 47.3%, sexM
    # 17005-17253 at 31.33%, fully contained. Ground truth is Plus.
    hits = [_hit("sexP", 16939, 17547, 47.3), _hit("sexM", 17005, 17253, 31.33)]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    by_gene = _by_gene(resolved)
    assert by_gene["sexM"].superseded_by == "sexP"
    assert by_gene["sexP"].superseded_by is None
    assert len(events) == 1


def test_the_superseded_hit_is_kept_not_dropped():
    # The ambiguity must stay visible in the report, and the observations are
    # what the 50% overlap threshold will later be recalibrated against.
    hits = [_hit("sexP", 16939, 17547, 47.3), _hit("sexM", 17005, 17253, 31.33)]
    resolved, _ = resolve_idiomorph_overlaps(hits, FAM)
    assert len(resolved) == 2
    assert {h.gene_name for h in resolved} == {"sexP", "sexM"}


def test_a_thin_margin_still_resolves_and_reports_its_margin():
    # Cunninghamella bertholletiae NRRL 1376: sexP 28.26 vs sexM 30.60, the
    # narrowest call in the ground-truth set, and correct (truth is Minus).
    # The call is still made; the margin is reported so the tier can be capped.
    hits = [_hit("sexP", 100, 400, 28.26), _hit("sexM", 120, 380, 30.60)]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert _by_gene(resolved)["sexP"].superseded_by == "sexM"
    assert events[0].winner == "sexM"
    assert events[0].loser == "sexP"
    assert abs(events[0].margin - 2.34) < 1e-9


def test_non_overlapping_mutually_exclusive_genes_are_both_kept():
    # Two DISTINCT proteins, one per idiomorph: a genuine both-idiomorphs
    # locus. Collapsing this would make that biology undetectable.
    hits = [_hit("sexP", 100, 400, 47.3), _hit("sexM", 5000, 5300, 44.0)]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert all(h.superseded_by is None for h in resolved)
    assert events == []


def test_overlap_below_the_threshold_does_not_resolve():
    # 300 bp hits sharing 100 bp: 33% of the shorter, under the 50% bar.
    hits = [_hit("sexP", 100, 399, 47.3), _hit("sexM", 300, 599, 31.3)]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert all(h.superseded_by is None for h in resolved)
    assert events == []


def test_overlap_is_measured_against_the_shorter_hit():
    # A short hit fully inside a long one overlaps it 100%, even though the
    # shared span is a small fraction of the LONGER hit. Measuring against the
    # longer hit would leave the real sexM/sexP case unresolved: in Absidia the
    # 249 bp sexM hit is only 41% of the 609 bp sexP hit.
    hits = [_hit("sexP", 1000, 1608, 47.3), _hit("sexM", 1066, 1314, 31.33)]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert _by_gene(resolved)["sexM"].superseded_by == "sexP"
    assert abs(events[0].overlap_fraction - 1.0) < 1e-9


def test_hits_on_different_contigs_never_resolve():
    hits = [_hit("sexP", 100, 400, 47.3), _hit("sexM", 100, 400, 31.3, contig="c2")]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert all(h.superseded_by is None for h in resolved)
    assert events == []


def test_genes_that_are_not_mutually_exclusive_are_never_resolved():
    # tptA has no present_in_idiomorphs, so it applies to every idiomorph.
    # A flanking gene overlapping a core gene is intergenic erosion -- normal
    # at MAT loci -- not an idiomorph ambiguity.
    hits = [
        _hit("tptA", 100, 400, 75.9, role="flanking_conserved"),
        _hit("sexP", 120, 380, 47.3),
    ]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert all(h.superseded_by is None for h in resolved)
    assert events == []


def test_two_hits_for_the_SAME_gene_are_never_resolved_against_each_other():
    # Gene duplication and multi-allele co-occurrence are normal at MAT loci.
    # Two sexP hits are two copies, not an idiomorph conflict.
    hits = [_hit("sexP", 100, 400, 47.3), _hit("sexP", 120, 380, 44.0)]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert all(h.superseded_by is None for h in resolved)
    assert events == []


def test_an_exact_identity_tie_is_left_unresolved():
    # With nothing to choose between them, picking one would be arbitrary and
    # the choice would flip on floating-point noise. Leaving both standing
    # reports the locus as undetermined, which is the honest answer.
    hits = [_hit("sexP", 100, 400, 35.0), _hit("sexM", 120, 380, 35.0)]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert all(h.superseded_by is None for h in resolved)
    assert events == []


def test_the_event_records_what_a_recalibration_would_need():
    hits = [
        _hit("sexP", 1000, 1608, 47.3, coverage=23.6),
        _hit("sexM", 1066, 1314, 31.33, coverage=None),
    ]
    _, events = resolve_idiomorph_overlaps(hits, FAM)
    event = events[0]
    assert event.contig == "c1"
    assert event.winner == "sexP"
    assert event.loser == "sexM"
    assert event.winner_identity == 47.3
    assert event.loser_identity == 31.33
    assert event.winner_coverage == 23.6
    assert event.loser_coverage is None


def test_a_custom_overlap_threshold_is_honoured():
    # 300 bp hits sharing 100 bp = 33% of the shorter. Under the 0.5 default
    # this does not resolve; at 0.3 it does. The threshold is provisional and
    # expected to be raised once the diagnostics corpus supports a value.
    hits = [_hit("sexP", 100, 399, 47.3), _hit("sexM", 300, 599, 31.3)]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM, min_overlap_fraction=0.3)
    assert _by_gene(resolved)["sexM"].superseded_by == "sexP"
    assert len(events) == 1
