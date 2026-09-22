"""`partial_locus` -- a call that cleared the admission bar only by tying it,
and every call the relaxed pass produces.

Curator's rulings, 2026-09-21:

* the class fires at `fraction_found == ambiguity_floor`, reusing the floor
  already in force rather than inventing a completeness constant. PROVISIONAL:
  to be revisited once there are examples to look at.
* it caps confidence at `medium`. `high` is the report's assertion "this is a
  MAT locus", and a call that only tied the floor must not make it.
* a relaxed-pass call NEVER emits `mat_locus`. Measured reason: on 283 BFD
  genomes the relaxed pass produced 352 loci at median identity 38.1% against
  strict's 77.6%, and one genome (Rhizomucor miehei CAU432) reported 13 loci,
  9 of them `mat_locus`, whose idiomorph calls contradicted each other
  (3 Plus / 4 Minus / 6 undetermined). `detection_pass` still distinguishes the
  two routes, so nothing is lost.
* the class displaces ONLY `mat_locus`. `idiomorph_gene_only`,
  `flanking_gene_only` and `homothallic_candidate` are statements about what
  genes are present, not weaker locus claims, so a fraction tie must not
  overwrite them.
"""
import pytest

from MATPredict.detect.idiomorph import (
    LOCUS_CLASS_FLANKING_ONLY,
    LOCUS_CLASS_HOMOTHALLIC,
    LOCUS_CLASS_IDIOMORPH_ONLY,
    LOCUS_CLASS_MAT,
    LOCUS_CLASS_PARTIAL,
    apply_partial_locus,
)
from MATPredict.detect.tiering import cap_at_medium


def test_a_mat_locus_that_ties_the_floor_becomes_partial():
    assert apply_partial_locus(
        LOCUS_CLASS_MAT, fraction_found=0.5, ambiguity_floor=0.5, relaxed=False
    ) == LOCUS_CLASS_PARTIAL


def test_a_mat_locus_above_the_floor_is_untouched():
    """The validated Mucoromycota calls sit at 0.833 (Plus) and 0.800 (Minus)
    and must stay `mat_locus`."""
    for fraction in (0.667, 0.800, 0.833, 1.0):
        assert apply_partial_locus(
            LOCUS_CLASS_MAT, fraction_found=fraction, ambiguity_floor=0.5, relaxed=False
        ) == LOCUS_CLASS_MAT


def test_the_rule_follows_a_moved_floor():
    """The class is defined against the floor in force, not against 0.5."""
    assert apply_partial_locus(
        LOCUS_CLASS_MAT, fraction_found=0.6, ambiguity_floor=0.6, relaxed=False
    ) == LOCUS_CLASS_PARTIAL
    assert apply_partial_locus(
        LOCUS_CLASS_MAT, fraction_found=0.5, ambiguity_floor=0.6, relaxed=False
    ) == LOCUS_CLASS_MAT


def test_every_relaxed_mat_locus_becomes_partial_whatever_its_fraction():
    for fraction in (0.333, 0.5, 0.667, 1.0):
        assert apply_partial_locus(
            LOCUS_CLASS_MAT, fraction_found=fraction, ambiguity_floor=0.5, relaxed=True
        ) == LOCUS_CLASS_PARTIAL


@pytest.mark.parametrize(
    "other",
    [LOCUS_CLASS_IDIOMORPH_ONLY, LOCUS_CLASS_FLANKING_ONLY, LOCUS_CLASS_HOMOTHALLIC],
)
@pytest.mark.parametrize("relaxed", [False, True])
def test_other_classes_are_never_displaced(other, relaxed):
    """Composition classes say which genes are present. A fraction tie, or the
    relaxed route, says how strong the call is -- a different axis."""
    assert apply_partial_locus(
        other, fraction_found=0.5, ambiguity_floor=0.5, relaxed=relaxed
    ) == other


def test_cap_at_medium_only_lowers():
    assert cap_at_medium("high") == "medium"
    assert cap_at_medium("medium") == "medium"
    assert cap_at_medium("low") == "low"


def test_cap_at_medium_closes_the_measured_hole():
    """All 5 families that can be `mat_locus` at exactly 0.500 carry a
    `flanking_conserved` gene, so `tptA + sexP + sexM` = 3/6 = 0.500 in
    Mucoromycota MAT satisfies core_found AND flanking_found and reaches
    `high`. The 8 spurious calls measured on the ground-truth set dodged this
    only because `glrA` is flanking_variable."""
    assert cap_at_medium("high") == "medium"


# --- the chain, against the REAL curated roster -----------------------------
# The unit tests above exercise each rule alone. This one runs scoring ->
# classification -> tiering together on the live db/Mucoromycota/order.yml, so
# a wiring mistake between them cannot pass unnoticed.

def _real_mucoro_family():
    from pathlib import Path
    from MATPredict.detect.family_registry import load_all_families
    return next(
        f for f in load_all_families(Path("db"))
        if f.key.phylum == "Mucoromycota" and f.key.locus_name == "MAT"
    )


def _cluster_of(family, spec):
    from MATPredict.detect.clustering import GeneCluster
    from MATPredict.detect.search import SearchHit
    roles = {g["name"]: g["role"] for g in family.genes}
    hits = [
        SearchHit(family.key, name, roles[name], "c1", start, start + 500, "+",
                  80.0, "rec1", method)
        for name, start, method in spec
    ]
    return GeneCluster("c1", min(h.start for h in hits), max(h.end for h in hits), hits)


def test_tptA_sexP_sexM_is_partial_and_medium_on_the_real_roster():
    """The hole identified on 2026-09-21: `tptA + sexP + sexM` is 3 of the 6
    non-optional Mucoromycota MAT genes = exactly 0.500, and because `tptA` is
    `flanking_conserved` it satisfies core_found AND flanking_found -- so
    before this change it scored `high`. The 8 spurious calls measured on the
    ground-truth set escaped `high` only because `glrA` is `flanking_variable`.

    `sexM` arrives by tblastn here, which is what actually happens: the
    proteome carries the gene that is really present and the genome-wide
    rescue lands the OTHER idiomorph's reference on it. Two proteome hits
    would instead be the documented homothallic architecture -- covered below.
    """
    from MATPredict.detect.idiomorph import classify_locus
    from MATPredict.detect.scoring import score_cluster
    from MATPredict.detect.tiering import assign_tier

    family = _real_mucoro_family()
    cluster = _cluster_of(family, [
        ("tptA", 1000, "diamond_proteome"),
        ("sexP", 3000, "diamond_proteome"),
        ("sexM", 5000, "tblastn_genome"),
    ])
    score = next(s for s in score_cluster(cluster, [family]) if s.family_key == family.key)
    assert score.fraction_found == 0.5, "roster changed; this test's premise is gone"

    tier = assign_tier(score, family, cluster, any_gene_unpolished=False, fragmented=False)
    assert tier == "high", "the hole this rule closes no longer exists"

    assert classify_locus(cluster, family) == LOCUS_CLASS_MAT
    assert apply_partial_locus(
        LOCUS_CLASS_MAT, score.fraction_found, 0.5, relaxed=False
    ) == LOCUS_CLASS_PARTIAL
    assert cap_at_medium(tier) == "medium"


def test_a_homothallic_candidate_at_the_floor_keeps_its_class_but_not_high():
    """Curator's ruling: the class is a statement about which genes are
    present and must survive; `high` is a statement about certainty and must
    not. Both sexP and sexM from the proteome on one contig is the documented
    homothallic architecture, so the class is correct -- but it rests on
    exactly half the roster."""
    from MATPredict.detect.idiomorph import classify_locus, is_partial_strength
    from MATPredict.detect.scoring import score_cluster

    family = _real_mucoro_family()
    cluster = _cluster_of(family, [
        ("tptA", 1000, "diamond_proteome"),
        ("sexP", 3000, "diamond_proteome"),
        ("sexM", 5000, "diamond_proteome"),
    ])
    score = next(s for s in score_cluster(cluster, [family]) if s.family_key == family.key)
    assert score.fraction_found == 0.5
    assert classify_locus(cluster, family) == LOCUS_CLASS_HOMOTHALLIC

    # class survives ...
    assert apply_partial_locus(
        LOCUS_CLASS_HOMOTHALLIC, score.fraction_found, 0.5, relaxed=False
    ) == LOCUS_CLASS_HOMOTHALLIC
    # ... but the cap still applies, because it keys on strength, not class.
    assert is_partial_strength(score.fraction_found, 0.5, relaxed=False) is True
