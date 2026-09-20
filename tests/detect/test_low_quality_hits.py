# tests/detect/test_low_quality_hits.py
"""Alignments too short to mean anything are dropped before they are counted.

Found in the 44-genus sweep on Syzygites sp. MES_3091, which reported as
separate "loci":

    scaffold_100  sexM 26,628-26,675 (48 bp)  sexP 26,643-26,693 (51 bp)
    scaffold_190  sexM  7,839-7,865 (27 bp)

A 27 bp alignment is nine codons. At that length a 77.8% identity is noise, not
evidence, and two such fragments 32 bp apart were being reported as a locus
with two genes. The curator ruled on 2026-09-20 that these are low-quality
candidates to be filtered rather than classified.

The floor is on the ALIGNMENT, not on the gene: a real MAT gene can be short
(the pheromone precursors this project exists to find are ~60-80 aa), but a
real hit to one covers a meaningful part of its reference rather than 9 codons.
"""
from __future__ import annotations

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.search import MIN_HIT_LENGTH_BP, SearchHit, drop_low_quality_hits

FAM = Family(
    FamilyKey("Mucoromycota", "MAT"), "enum", ["Plus", "Minus"], None,
    [
        {"name": "sexP", "role": "core_MAT", "present_in_idiomorphs": ["Plus"]},
        {"name": "sexM", "role": "core_MAT", "present_in_idiomorphs": ["Minus"]},
    ],
    [4827],
)


def _hit(gene, start, end, identity=50.0, method="tblastn_genome"):
    return SearchHit(
        FAM.key, gene, "core_MAT", "c1", start, end, "+", identity, "rec1", method
    )


def test_a_twenty_seven_base_alignment_is_dropped():
    # The real Syzygites scaffold_190 sexM hit.
    kept = drop_low_quality_hits([_hit("sexM", 7839, 7865, identity=77.778)])
    assert kept == []


def test_the_fifty_base_syzygites_pair_is_dropped_entirely():
    kept = drop_low_quality_hits([
        _hit("sexM", 26628, 26675, identity=50.0),
        _hit("sexP", 26643, 26693, identity=52.941),
    ])
    assert kept == []


def test_a_hit_at_the_floor_is_kept():
    kept = drop_low_quality_hits([_hit("sexM", 1000, 1000 + MIN_HIT_LENGTH_BP - 1)])
    assert len(kept) == 1


def test_a_full_length_hit_is_kept_however_low_its_identity():
    # Identity is NOT the filter. Minus-strain identities in this project's own
    # ground-truth set run 25.9-43.5%, so an identity floor would destroy Minus
    # detection -- the reason min_identity is still None.
    kept = drop_low_quality_hits([_hit("sexM", 1000, 1600, identity=25.9)])
    assert len(kept) == 1


def test_a_short_but_real_pheromone_gene_hit_is_kept():
    # 60 aa = 180 bp, the short-ORF floor this project already uses. The length
    # bar must sit well below that or it would drop the small MAT genes the
    # pipeline exists to rescue.
    assert MIN_HIT_LENGTH_BP < 180
    kept = drop_low_quality_hits([_hit("sexM", 1000, 1179)])
    assert len(kept) == 1
