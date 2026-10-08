"""`supported_span`: a locus's span rebuilt from its own modelled genes plus its own hits strong enough to trust.

Report only. The cluster span chains every hit of every family at BLAST's default e-value, so one weak hit can stretch
a locus without adding a gene (T48-F: a called HD locus of 12.3 kb reported as 35.2 kb because of one hit with
e-value 4.1; analysis/2026-10-08_cinerea-b43-trace.md). Trimming to the modelled genes would also drop real,
unmodelled receptor genes (6 of 11 loci with 10 kb or more beyond their genes), so this keeps own-family hits at or
above a bitscore floor. A bitscore, not an e-value: an e-value depends on genome size, and the flank-carried rule
already moved to a bitscore floor for that reason (curator, 2026-09-27).
"""
from __future__ import annotations

DEFAULT_SUPPORTED_MIN_BITSCORE = 39.0   # the flank-carried floor (family_registry.DEFAULT_FLANK_CARRIED_MIN_BITSCORE)


def supported_span(clusters, family_key, contig, evidence, span_start: int, span_end: int,
                   min_bitscore: float = DEFAULT_SUPPORTED_MIN_BITSCORE) -> dict | None:
    """Extent on `contig` of (a) the locus's modelled genes and (b) its own-family hits with bitscore >= `min_bitscore`.

    A hit with no bitscore (a polished model, a CAAX scan hit) is not counted as support by itself; polished models
    count through `evidence`. A hit that lost an idiomorph resolution is skipped. `beyond_supported_bp` is how much of
    the reported span (`span_start`-`span_end`) lies outside the supported extent. None when nothing lies on `contig`."""
    points = [(g.start, g.end) for g in evidence if g.contig == contig]
    for cluster in clusters:
        for h in cluster.hits:
            if (h.family_key == family_key and h.contig == contig and h.bitscore is not None
                    and h.bitscore >= min_bitscore and getattr(h, "superseded_by", None) is None):
                points.append((h.start, h.end))
    if not points:
        return None
    start, end = min(p[0] for p in points), max(p[1] for p in points)
    return {"start": start, "end": end, "min_bitscore": min_bitscore,
            "beyond_supported_bp": (span_end - span_start + 1) - (end - start + 1)}
