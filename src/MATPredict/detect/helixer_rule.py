"""The flank-only helixer rule: which gene model to report when homology and helixer-lite disagree.

Decided from the panel metric (docs/notes/2026-10-15_session-handoff.md, "Metric results"): helixer-lite's
ends agree with the NCBI annotation more often than homology's on flanking genes, and less often on core MAT
genes, where it also undercalls. So:

* A **flanking** gene (role ``flanking_*``) takes helixer-lite's ends when a helixer model matches it (same
  contig and strand, overlap at least half the longer of the two) and that model is not known to be an
  incomplete ORF.
* Every other gene (``core_MAT`` and the rest) keeps its homology model. A missing helixer model is no
  evidence against homology; silence is not contradiction.
* The rule never removes or adds a gene and never changes a call. It returns the refined model beside the
  original, with the action and the reason, so the change is recorded and the original stays in the report.

``complete_orf`` None means unknown (an older helixer table without the column) and does not block the rule;
only an explicit False does.
"""
from __future__ import annotations

from dataclasses import dataclass

from MATPredict.detect.tracks import Track, TrackGene

MIN_RECIPROCAL = 0.5

KEPT = "kept"
REFINED = "refined"


@dataclass(frozen=True)
class Refinement:
    gene: str
    role: str
    contig: str
    strand: str
    action: str  # KEPT or REFINED
    reason: str
    start: int  # the model to report
    end: int
    original_start: int
    original_end: int
    helixer_gene: str | None = None
    delta_5p: int = 0  # bp the reported end moved from the homology end, on the gene's strand; + is outward
    delta_3p: int = 0


def is_flanking(role: str | None) -> bool:
    return bool(role) and role.startswith("flanking")


def best_match(contig: str, start: int, end: int, strand: str, models) -> TrackGene | None:
    """The helixer model with the largest overlap that passes the reciprocal 50% rule, or None."""
    best, best_ov = None, 0
    for m in models:
        if m.contig != contig or m.strand != strand:
            continue
        ov = min(end, m.end) - max(start, m.start) + 1
        if ov <= 0:
            continue
        longer = max(end - start + 1, m.end - m.start + 1)
        if ov >= MIN_RECIPROCAL * longer and ov > best_ov:
            best, best_ov = m, ov
    return best


def _deltas(strand: str, start: int, end: int, new_start: int, new_end: int) -> tuple[int, int]:
    """(5' delta, 3' delta) in bp on the gene's strand; positive means the model grew at that end."""
    if strand == "+":
        return start - new_start, new_end - end
    return new_end - end, start - new_start


def refine_gene(gene, helixer: Track) -> Refinement:
    """Apply the rule to one homology gene. `gene` needs gene_name, role, contig, start, end and strand."""
    def kept(reason: str, hx: TrackGene | None = None) -> Refinement:
        return Refinement(gene.gene_name, gene.role, gene.contig, gene.strand, KEPT, reason,
                          gene.start, gene.end, gene.start, gene.end, hx.name if hx else None)

    if not is_flanking(gene.role):
        return kept("homology_stratum")
    hx = best_match(gene.contig, gene.start, gene.end, gene.strand, helixer.genes)
    if hx is None:
        return kept("no_helixer_model")
    if hx.complete_orf is False:
        return kept("helixer_orf_incomplete", hx)
    if (hx.start, hx.end) == (gene.start, gene.end):
        return kept("ends_agree", hx)
    d5, d3 = _deltas(gene.strand, gene.start, gene.end, hx.start, hx.end)
    return Refinement(gene.gene_name, gene.role, gene.contig, gene.strand, REFINED, "helixer_ends",
                      hx.start, hx.end, gene.start, gene.end, hx.name, d5, d3)


def refine_genes(genes, helixer: Track) -> list[Refinement]:
    """The rule over every present gene of a detection; one Refinement per gene, in input order."""
    return [refine_gene(g, helixer) for g in genes]
