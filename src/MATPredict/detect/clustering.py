"""Group search hits into spatial clusters by max intergenic gap (strand- and family-agnostic)."""
from __future__ import annotations

from dataclasses import dataclass

from MATPredict.detect.search import SearchHit


@dataclass(frozen=True)
class GeneCluster:
    contig: str
    start: int
    end: int
    hits: list[SearchHit]


def cluster_hits(hits: list[SearchHit], max_gap: int = 25_000) -> list[GeneCluster]:
    """Sort hits per contig by start, then split into clusters wherever the
    gap to the next hit's start exceeds max_gap. Never groups across contigs.
    Strand is deliberately not a grouping criterion (real curated loci mix
    strands within one locus)."""
    clusters: list[GeneCluster] = []
    by_contig: dict[str, list[SearchHit]] = {}
    for hit in hits:
        by_contig.setdefault(hit.contig, []).append(hit)

    for contig, contig_hits in by_contig.items():
        contig_hits.sort(key=lambda h: h.start)
        current: list[SearchHit] = [contig_hits[0]]
        current_end = contig_hits[0].end
        for hit in contig_hits[1:]:
            if hit.start - current_end > max_gap:
                clusters.append(GeneCluster(contig, current[0].start, current_end, current))
                current = [hit]
                current_end = hit.end
            else:
                current.append(hit)
                current_end = max(current_end, hit.end)
        clusters.append(GeneCluster(contig, current[0].start, current_end, current))
    return clusters
