"""Rescue polish is not run for the idiomorph the cluster demonstrably is not.

Curator's ruling, 2026-09-22 ("1 do this"), taken against a measured cost: one
rescue polish pair is exonerate 0.77 s + miniprot 0.12 s ~= 0.9 s, and a
genome-only Pezizomycotina run spends ~230 s in polish while only a handful of
those pairs yield evidence that is used. A cluster that already evidences
MAT1-1 cannot also hold MAT1-2-1, so rescuing for it is guaranteed-futile.

The narrowing is deliberately timid. It fires ONLY on exactly one evidenced
idiomorph; flanking-only clusters (zero evidenced -- the localize-by-flanks
case this project depends on) and genuine both-idiomorph loci (two or more)
keep the full rescue set.
"""
from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.idiomorph import evidenced_idiomorphs
from MATPredict.detect.pipeline import _rescue_genes_in_idiomorph_scope
from MATPredict.detect.search import SearchHit

KEY = FamilyKey("Ascomycota", "MAT")
FAM = Family(
    key=KEY, vocabulary_type="enum", idiomorph_values=["MAT1-1", "MAT1-2"],
    idiomorph_pattern=None,
    genes=[
        {"name": "MAT1-1-1", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-1"]},
        {"name": "MAT1-1-3", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-1"]},
        {"name": "MAT1-2-1", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-2"]},
        {"name": "APN2", "role": "flanking_conserved"},
    ],
    taxonomic_scope=[4890],
)
ALL_CORE = {"MAT1-1-1", "MAT1-1-3", "MAT1-2-1"}


def _hit(gene, role="core_MAT", superseded_by=None):
    return SearchHit(KEY, gene, role, "c1", 1000, 1300, "+", 80.0, "r",
                     "tblastn_genome", bitscore=200.0, superseded_by=superseded_by)


def _cluster(*hits):
    return GeneCluster(contig="c1", start=900, end=2000, hits=list(hits))


def test_the_other_idiomorphs_rescue_is_skipped():
    cluster = _cluster(_hit("MAT1-1-1"))
    keep, skip = _rescue_genes_in_idiomorph_scope(cluster, FAM, set(ALL_CORE))
    assert skip == {"MAT1-2-1"}
    assert keep == {"MAT1-1-1", "MAT1-1-3"}


def test_a_flanking_only_cluster_keeps_the_full_rescue_set():
    # The localize-by-flanks case: nothing informative yet, so nothing to
    # narrow by. Narrowing here would defeat the whole point of flank anchoring.
    cluster = _cluster(_hit("APN2", role="flanking_conserved"))
    keep, skip = _rescue_genes_in_idiomorph_scope(cluster, FAM, set(ALL_CORE))
    assert skip == set()
    assert keep == ALL_CORE


def test_a_both_idiomorph_cluster_keeps_the_full_rescue_set():
    cluster = _cluster(_hit("MAT1-1-1"), _hit("MAT1-2-1"))
    keep, skip = _rescue_genes_in_idiomorph_scope(cluster, FAM, set(ALL_CORE))
    assert skip == set()
    assert keep == ALL_CORE


def test_a_superseded_hit_does_not_narrow_the_rescue():
    # The losing half of a resolved cross-match is not evidence for its own
    # idiomorph, so a cluster whose only MAT1-2 hit was superseded is a
    # single-idiomorph cluster -- and one whose only MAT1-1 hit was superseded
    # evidences nothing at all and must not narrow.
    cluster = _cluster(_hit("MAT1-1-1", superseded_by="x"))
    keep, skip = _rescue_genes_in_idiomorph_scope(cluster, FAM, set(ALL_CORE))
    assert skip == set()
    assert keep == ALL_CORE


def test_an_idiomorph_agnostic_gene_is_never_skipped():
    cluster = _cluster(_hit("MAT1-1-1"))
    keep, skip = _rescue_genes_in_idiomorph_scope(
        cluster, FAM, {"MAT1-2-1", "APN2"}
    )
    assert skip == {"MAT1-2-1"}
    assert "APN2" in keep


def test_evidenced_idiomorphs_ignores_uninformative_genes():
    fam = Family(
        key=KEY, vocabulary_type="enum", idiomorph_values=["MAT1-1", "MAT1-2"],
        idiomorph_pattern=None,
        genes=[
            {"name": "MAT1-1-1", "role": "core_MAT",
             "present_in_idiomorphs": ["MAT1-1"], "idiomorph_informative": False},
        ],
        taxonomic_scope=[4890],
    )
    assert evidenced_idiomorphs(fam, [_hit("MAT1-1-1")]) == set()
