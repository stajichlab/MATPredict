# tests/detect/test_locus_class.py
"""Detected loci are classified, not filtered.

Curator's ruling, 2026-09-20, after the 44-genus sweep produced 81 calls whose
only genes were sexP and sexM: these are NOT junk to be thrown away.

* A lone sexM or lone sexP, with no flanking gene at all, is still valuable --
  it is training material for a per-idiomorph HMM, which is a search strategy
  this project intends to build. It belongs in its own category rather than
  being discarded or mixed in with confirmed loci.
* BOTH idiomorphs in one locus is the documented homothallic architecture.
  Schulz et al. 2016 and Idnurm 2011 describe Syzygites megalocarpus, a
  homothallic Mucorale, as encoding both HMG transcription factors, each
  flanked by its own glrA or rnhA. Collapsing or rejecting those would make
  homothallism undetectable.

So `locus_class` describes WHAT was found; `detection_pass` describes HOW it
was admitted. They are orthogonal on purpose -- a homothallic candidate can
come from either pass.
"""
from __future__ import annotations

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import (
    DEFAULT_MAX_HOMOTHALLIC_SEPARATION_BP,
    Family,
    FamilyKey,
)
from MATPredict.detect.idiomorph import classify_locus
from MATPredict.detect.search import SearchHit

FAM = Family(
    FamilyKey("Mucoromycota", "MAT"), "enum", ["Plus", "Minus"], None,
    [
        {"name": "tptA", "role": "flanking_conserved"},
        {"name": "sexP", "role": "core_MAT", "present_in_idiomorphs": ["Plus"]},
        {"name": "sexM", "role": "core_MAT", "present_in_idiomorphs": ["Minus"]},
        {"name": "rnhA", "role": "flanking_conserved"},
        {"name": "glrA", "role": "flanking_variable"},
    ],
    [4827],
)


def _hit(gene, start, end, role="core_MAT", contig="c1"):
    return SearchHit(
        FAM.key, gene, role, contig, start, end, "+", 60.0, "rec1",
        "diamond_proteome",
    )


def _cluster(*hits):
    return GeneCluster(hits[0].contig, min(h.start for h in hits),
                       max(h.end for h in hits), list(hits))


def test_a_core_gene_with_a_flank_is_an_ordinary_mat_locus():
    cluster = _cluster(
        _hit("sexP", 1000, 1600),
        _hit("rnhA", 2000, 5000, role="flanking_conserved"),
    )
    assert classify_locus(cluster, FAM) == "mat_locus"


def test_a_lone_idiomorph_gene_is_its_own_category_not_a_locus():
    # Kept deliberately: this is HMM training material for the idiomorph even
    # with no linked flanking gene.
    cluster = _cluster(_hit("sexM", 1000, 1600))
    assert classify_locus(cluster, FAM) == "idiomorph_gene_only"


def test_two_copies_of_one_idiomorph_gene_are_still_idiomorph_gene_only():
    cluster = _cluster(_hit("sexM", 1000, 1600), _hit("sexM", 8000, 8600))
    assert classify_locus(cluster, FAM) == "idiomorph_gene_only"


def test_both_idiomorphs_close_together_is_a_homothallic_candidate():
    # The Syzygites architecture: both HMG genes in one locus.
    cluster = _cluster(_hit("sexP", 1000, 1600), _hit("sexM", 4000, 4600))
    assert classify_locus(cluster, FAM) == "homothallic_candidate"


def test_both_idiomorphs_with_flanks_is_still_a_homothallic_candidate():
    # Flanks make the homothallic call STRONGER, not weaker -- Schulz describes
    # each HMG gene flanked by its own glrA or rnhA.
    cluster = _cluster(
        _hit("sexP", 1000, 1600),
        _hit("rnhA", 2000, 5000, role="flanking_conserved"),
        _hit("sexM", 6000, 6600),
        _hit("glrA", 7000, 9000, role="flanking_variable"),
    )
    assert classify_locus(cluster, FAM) == "homothallic_candidate"


def test_idiomorphs_too_far_apart_are_not_a_homothallic_candidate():
    # Beyond the separation bar they are two unrelated idiomorph-gene hits that
    # the clustering gap happened to group, not one homothallic locus.
    far = DEFAULT_MAX_HOMOTHALLIC_SEPARATION_BP + 5000
    cluster = _cluster(_hit("sexP", 1000, 1600), _hit("sexM", far, far + 600))
    assert classify_locus(cluster, FAM) == "idiomorph_gene_only"


def test_idiomorphs_on_different_contigs_are_not_a_homothallic_candidate():
    hits = [_hit("sexP", 1000, 1600), _hit("sexM", 1000, 1600, contig="c2")]
    cluster = GeneCluster("c1", 1000, 1600, hits)
    assert classify_locus(cluster, FAM) == "idiomorph_gene_only"


def test_a_superseded_hit_does_not_make_a_homothallic_candidate():
    # A resolved cross-hit is ONE gene seen twice, not two idiomorphs present.
    # Without this, every ordinary heterothallic locus whose sexM/sexP overlap
    # was collapsed would be mislabelled homothallic.
    import dataclasses

    cluster = _cluster(
        _hit("sexP", 1000, 1600),
        dataclasses.replace(_hit("sexM", 1050, 1550), superseded_by="sexP"),
        _hit("rnhA", 2000, 5000, role="flanking_conserved"),
    )
    assert classify_locus(cluster, FAM) == "mat_locus"


def test_the_separation_default_is_the_curators_provisional_value():
    assert DEFAULT_MAX_HOMOTHALLIC_SEPARATION_BP == 20_000


def test_the_separation_is_measured_between_each_genes_BEST_hit():
    # A cluster holds many tblastn HSPs of the same region -- 15 or more is
    # routine. Testing every pairwise combination means some sexP/sexM pair is
    # almost always within any threshold, so the class stops discriminating:
    # the 44-genus sweep produced 77 homothallic candidates across 29 genera,
    # including separations of 43,879 bp under a 20 kb bar, because a stray
    # HSP pair was close even though each gene's best hit was far apart.
    #
    # The separation must be measured between the hits the report actually
    # shows as that gene's evidence -- the best one per gene -- so the class
    # and the displayed coordinates cannot disagree.
    far = DEFAULT_MAX_HOMOTHALLIC_SEPARATION_BP + 20_000
    cluster = _cluster(
        # Each gene's BEST hit (proteome path) is far apart...
        SearchHit(FAM.key, "sexP", "core_MAT", "c1", 1000, 1600, "+", 60.0,
                  "rec1", "diamond_proteome"),
        SearchHit(FAM.key, "sexM", "core_MAT", "c1", far, far + 600, "+", 60.0,
                  "rec1", "diamond_proteome"),
        # ...but a weak stray HSP of sexM sits right next to sexP.
        SearchHit(FAM.key, "sexM", "core_MAT", "c1", 2000, 2200, "+", 22.0,
                  "rec1", "tblastn_genome"),
    )
    assert classify_locus(cluster, FAM) == "idiomorph_gene_only"


def test_best_hits_close_together_still_call_homothallic():
    cluster = _cluster(
        SearchHit(FAM.key, "sexP", "core_MAT", "c1", 1000, 1600, "+", 60.0,
                  "rec1", "diamond_proteome"),
        SearchHit(FAM.key, "sexM", "core_MAT", "c1", 4000, 4600, "+", 58.0,
                  "rec1", "diamond_proteome"),
        SearchHit(FAM.key, "sexM", "core_MAT", "c1", 90000, 90200, "+", 20.0,
                  "rec1", "tblastn_genome"),
    )
    assert classify_locus(cluster, FAM) == "homothallic_candidate"
