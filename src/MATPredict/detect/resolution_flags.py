"""`resolution_flags`: why a called locus may not be a resolved answer (report only; calls, scores and tiers are unchanged).

A locus can be called with confidence and still be a fragment: its core gene runs off the end of a contig, covers a small part of its
reference, or is only a few dozen residues long. These flags say so next to the call (curator ruling 2026-11-04, from the divergent
Ogataea MATog loci, #256) and never feed back into it. Seam: the detect report writer.

  contig_edge   a core gene lies within EDGE_BP of a contig end (needs the contig lengths)
  partial_core  a core gene covers less than PARTIAL_COVERAGE % of its reference protein
  short_core    a core gene's model is shorter than SHORT_CORE_AA residues
"""
from __future__ import annotations

EDGE_BP = 500
PARTIAL_COVERAGE = 50.0
SHORT_CORE_AA = 80
FLAGS = ("contig_edge", "partial_core", "short_core")


def _model_aa(e) -> int:
    """Residues in the gene model: the exons when listed, else the gene extent, divided by 3."""
    exons = getattr(e, "exons", None)
    n = sum(b - a + 1 for a, b in exons) if exons else e.end - e.start + 1
    return n // 3


def resolution_flags(gene_evidence, contig_lengths: dict | None = None) -> dict[str, list[str]]:
    """{flag: [gene names]} for the locus's own core_MAT genes; a flag with no gene is absent. The losing half of a resolved pair is not this
    locus's gene and is skipped."""
    out: dict[str, set[str]] = {}
    for e in gene_evidence or []:
        if e.role != "core_MAT" or getattr(e, "lost_to", None):
            continue
        length = (contig_lengths or {}).get(e.contig)
        if length is not None and min(e.start - 1, length - e.end) < EDGE_BP:
            out.setdefault("contig_edge", set()).add(e.gene_name)
        if e.coverage is not None and e.coverage < PARTIAL_COVERAGE:
            out.setdefault("partial_core", set()).add(e.gene_name)
        if _model_aa(e) < SHORT_CORE_AA:
            out.setdefault("short_core", set()).add(e.gene_name)
    return {k: sorted(out[k]) for k in FLAGS if k in out}
