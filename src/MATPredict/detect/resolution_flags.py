"""`resolution_flags`: why a called locus may not be a resolved answer (report only; calls, scores and tiers are unchanged).

A locus can be called with confidence and still be a fragment: its core gene runs off the end of a contig, covers a small part of its
reference, or is only a few dozen residues long. These flags say so next to the call (curator ruling 2026-11-04, from the divergent
Ogataea MATog loci, #256) and never feed back into it. Seam: the detect report writer.

  contig_edge   a core gene lies within EDGE_BP of a contig end (needs the contig lengths)
  partial_core  a core gene's model covers less than PARTIAL_COVERAGE % of its reference protein AND matches it at less than
                PARTIAL_IDENTITY % identity: a weak, partial model. A short exon piece at high identity is an intact intron-containing
                gene modelled in pieces (MATsc MATALPHA1, 40 aa at 97% identity) and is not flagged. A chromosome-level run
                (80 Saccharomycetaceae/Aspergillaceae/Nectriaceae genomes) showed coverage or length alone flagging 60 of 98 resolved loci.
"""
from __future__ import annotations

EDGE_BP = 500
PARTIAL_COVERAGE = 50.0
PARTIAL_IDENTITY = 90.0
FLAGS = ("contig_edge", "partial_core")


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
        if e.coverage is not None and e.coverage < PARTIAL_COVERAGE and e.identity is not None and e.identity < PARTIAL_IDENTITY:
            out.setdefault("partial_core", set()).add(e.gene_name)
    return {k: sorted(out[k]) for k in FLAGS if k in out}
