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


def test_the_proteome_found_gene_wins_even_at_lower_identity():
    # The discriminator is WHICH SEARCH FOUND IT, not identity. The fast path
    # matches the annotated proteome with diamond, so a gene the annotation
    # actually predicted is found there. The localization rescue then runs
    # tblastn genome-wide for core genes MISSING from the cluster, which lands
    # the other idiomorph's gene on the same spot precisely BECAUSE it is not
    # really there. The proteome hit is therefore the real gene.
    #
    # Measured: Cunninghamella bertholletiae NRRL 1376, a Minus genome, after
    # three published references were added. sexP rose to 34.146 (tblastn,
    # against the new Mooraboolomyces Plus reference) and overtook sexM at
    # 30.60 (diamond, against the annotated gene), flipping a correct call.
    # Identity scores 19/23 on that reference set; this rule scores 22/23, and
    # neither is worse on the smaller set.
    hits = [
        SearchHit(FAM.key, "sexP", "core_MAT", "c1", 100, 400, "+", 34.146,
                  "rec_P", "tblastn_genome"),
        SearchHit(FAM.key, "sexM", "core_MAT", "c1", 120, 380, "+", 30.60,
                  "rec_M", "diamond_proteome", coverage=94.8),
    ]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert _by_gene(resolved)["sexP"].superseded_by == "sexM"
    assert events[0].winner == "sexM"


def test_identity_decides_when_both_came_from_the_proteome():
    hits = [_hit("sexP", 100, 400, 47.3), _hit("sexM", 120, 380, 31.33)]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert _by_gene(resolved)["sexM"].superseded_by == "sexP"
    assert events[0].winner == "sexP"


def test_identity_decides_when_neither_came_from_the_proteome():
    # Both rescued by tblastn: the annotation predicted neither gene, so the
    # path carries no information and identity is all that is left. This is
    # the one ground-truth genome the rule cannot call (Cunninghamella
    # polymorpha NRRL 1395), and it is also the annotation-gap case this
    # project already knows about for small MAT genes.
    hits = [
        SearchHit(FAM.key, "sexP", "core_MAT", "c1", 100, 400, "+", 34.884,
                  "rec_P", "tblastn_genome"),
        SearchHit(FAM.key, "sexM", "core_MAT", "c1", 120, 380, "+", 33.929,
                  "rec_M", "tblastn_genome"),
    ]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert _by_gene(resolved)["sexM"].superseded_by == "sexP"


def test_a_polished_hit_still_counts_as_proteome_evidence():
    # A gene first found by diamond and then refined by exonerate/miniprot
    # carries the refiner's method, not diamond's. It must not lose its
    # proteome standing just because it was polished.
    hits = [
        SearchHit(FAM.key, "sexP", "core_MAT", "c1", 100, 400, "+", 34.0,
                  "rec_P", "tblastn_genome"),
        SearchHit(FAM.key, "sexM", "core_MAT", "c1", 120, 380, "+", 30.0,
                  "rec_M", "exonerate_refine"),
    ]
    resolved, _ = resolve_idiomorph_overlaps(hits, FAM)
    assert _by_gene(resolved)["sexP"].superseded_by == "sexM"


def test_the_default_overlap_bar_is_the_calibrated_value():
    # Briefly raised to 0.8 on the 23-genome corpus, then RETURNED to 0.5 when
    # the 44-genus sweep showed what the higher bar costs: a partial overlap of
    # ONE HMG region escapes collapse, so one gene is reported as two. Twelve of
    # 75 homothallic candidates had separations of zero or less -- overlapping,
    # and still reported as a pair. The 23 real loci overlap at 0.9242 or more,
    # so 0.5 loses nothing there either. Pinned so a future change is a
    # deliberate recalibration, not a drift.
    from MATPredict.detect.idiomorph import DEFAULT_MIN_OVERLAP_FRACTION

    assert DEFAULT_MIN_OVERLAP_FRACTION == 0.5


def test_a_partial_overlap_of_one_hmg_region_is_collapsed():
    # 1000 bp hits sharing 700 bp of the shorter = 0.70. This is the case the
    # 0.8 bar let through: two references aligning to overlapping but not
    # congruent spans of ONE gene, reported as two genes and mislabelled a
    # homothallic pair. At 0.5 it collapses.
    hits = [_hit("sexP", 1000, 1999, 47.3), _hit("sexM", 1300, 2299, 31.3)]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert _by_gene(resolved)["sexM"].superseded_by == "sexP"
    assert len(events) == 1


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


def test_several_references_per_gene_produce_one_verdict_not_a_pairwise_mess():
    # The curated Mucoromycota database holds 3 sexP and 3 sexM proteins, so
    # ONE real locus gene draws several hits of each name. Comparing every
    # sexP against every sexM produces contradictory verdicts -- on the real
    # Absidia cuneospora locus, 9 events of which one had sexM beating sexP --
    # and a spurious margin of 0.825 that wrongly capped the tier, when the
    # honest comparison of best against best is 47.3 vs 31.325.
    hits = [
        _hit("sexP", 16939, 17547, 47.3),
        _hit("sexP", 16939, 17547, 35.6),
        _hit("sexP", 16939, 17547, 30.5),
        _hit("sexM", 17005, 17253, 31.325),
        _hit("sexM", 17005, 17253, 29.63),
        _hit("sexM", 17005, 17253, 24.638),
    ]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert len(events) == 1
    assert events[0].winner == "sexP"
    assert events[0].winner_identity == 47.3
    assert events[0].loser_identity == 31.325
    assert abs(events[0].margin - 15.975) < 1e-9
    # EVERY losing hit is superseded, and no winning hit is.
    assert all(h.superseded_by == "sexP" for h in resolved if h.gene_name == "sexM")
    assert all(h.superseded_by is None for h in resolved if h.gene_name == "sexP")


def test_two_independent_loci_each_get_their_own_verdict():
    # Gene duplication is normal at MAT loci, so two non-overlapping groups
    # must be resolved separately rather than pooled into one comparison.
    hits = [
        _hit("sexP", 1000, 1600, 47.3), _hit("sexM", 1050, 1300, 31.3),
        _hit("sexP", 9000, 9600, 20.0), _hit("sexM", 9050, 9300, 44.0),
    ]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM)
    assert {(e.winner, e.loser) for e in events} == {("sexP", "sexM"), ("sexM", "sexP")}
    by_start = {(h.gene_name, h.start): h.superseded_by for h in resolved}
    assert by_start[("sexM", 1050)] == "sexP"   # first locus: sexP wins
    assert by_start[("sexP", 1000)] is None
    assert by_start[("sexP", 9000)] == "sexM"   # second locus: sexM wins
    assert by_start[("sexM", 9050)] is None


def test_a_custom_overlap_threshold_is_honoured():
    # 300 bp hits sharing 100 bp = 33% of the shorter. Under the 0.5 default
    # this does not resolve; at 0.3 it does. The threshold is provisional and
    # expected to be raised once the diagnostics corpus supports a value.
    hits = [_hit("sexP", 100, 399, 47.3), _hit("sexM", 300, 599, 31.3)]
    resolved, events = resolve_idiomorph_overlaps(hits, FAM, min_overlap_fraction=0.3)
    assert _by_gene(resolved)["sexM"].superseded_by == "sexP"
    assert len(events) == 1
