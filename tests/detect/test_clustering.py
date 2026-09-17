from __future__ import annotations

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.search import SearchHit
from MATPredict.detect.clustering import cluster_hits


def _hit(contig, start, end, strand="+"):
    return SearchHit(FamilyKey("P", "L"), "geneX", "core_MAT", contig, start, end, strand, 90.0, "rec1", "diamond_proteome")


def test_cluster_hits_groups_within_max_gap():
    hits = [_hit("c1", 100, 200), _hit("c1", 5000, 5100)]
    clusters = cluster_hits(hits, max_gap=10_000)
    assert len(clusters) == 1
    assert clusters[0].start == 100 and clusters[0].end == 5100


def test_cluster_hits_splits_beyond_max_gap():
    hits = [_hit("c1", 100, 200), _hit("c1", 50_000, 50_100)]
    clusters = cluster_hits(hits, max_gap=10_000)
    assert len(clusters) == 2


def test_cluster_hits_ignores_strand_for_grouping():
    hits = [_hit("c1", 100, 200, strand="+"), _hit("c1", 300, 400, strand="-")]
    clusters = cluster_hits(hits, max_gap=10_000)
    assert len(clusters) == 1


def test_cluster_hits_never_groups_across_contigs():
    hits = [_hit("c1", 100, 200), _hit("c2", 150, 250)]
    clusters = cluster_hits(hits, max_gap=10_000)
    assert len(clusters) == 2
