"""The flank-only helixer rule: which gene model to report when homology and helixer-lite disagree.

Decided from the panel metric (docs/notes/2026-10-15_session-handoff.md, "Metric results"): helixer-lite's
ends agree with the NCBI annotation more often than homology's on flanking genes, and less often on core MAT
genes, where it also undercalls. So:

* A **flanking** gene (role ``flanking_*``) is **extended** to helixer-lite's ends when a helixer model matches
  it (same contig and strand, overlap at least half the longer of the two) and that model is not known to be
  an incomplete ORF. Each end is taken on its own and only ever moves outward: where helixer stops short of
  the homology end, the homology end stays. Shortening is not allowed until a larger gene set (Aspergillus,
  Neurospora) has been tested; a first version that also shortened lost agreement on Aspergillus.
* Every other gene (``core_MAT`` and the rest) keeps its homology model. A missing helixer model is no
  evidence against homology; silence is not contradiction.
* The rule never removes or adds a gene and never changes a call. It returns the refined model beside the
  original, with the action and the reason, so the change is recorded and the original stays in the report.

**Query check.** Given the curated query protein, the homology protein and a way to translate a helixer model,
the extension is kept only when helixer's protein is better explained by the query than the homology protein
is: its local-alignment score against the query (BLOSUM62) must exceed the homology protein's by at least
`MIN_SCORE_GAIN`. An extension whose extra residues do not align to the query, or that replaces aligned
residues with unaligned ones, is not evidence of a better model and the homology model stays. Without the
three inputs the rule runs unverified (`Refinement.verified` False).

``complete_orf`` None means unknown (an older helixer table without the column) and does not block the rule;
only an explicit False does.
"""
from __future__ import annotations

from dataclasses import dataclass

from MATPredict.detect.tracks import Track, TrackGene

MIN_RECIPROCAL = 0.5
#: Raw BLOSUM62 local-alignment score helixer's protein must add over the homology protein (about three
#: identical residues). Fixed before the panel was scored with it.
MIN_SCORE_GAIN = 10.0

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
    verified: bool = False  # True when the query check ran
    score_gain: float | None = None  # helixer protein's query score minus the homology protein's


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


_ALIGNER = None


def query_score(query: str, protein: str) -> float:
    """Best local-alignment score (BLOSUM62, gap open -11, extend -1) of `protein` against `query`."""
    global _ALIGNER
    from Bio import Align
    from Bio.Align import substitution_matrices

    if _ALIGNER is None:
        a = Align.PairwiseAligner()
        a.mode = "local"
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score, a.extend_gap_score = -11, -1
        _ALIGNER = a
    q, p = query.replace("*", "").upper(), protein.replace("*", "").upper()
    return float(_ALIGNER.score(q, p)) if q and p else 0.0


def refine_gene(gene, helixer: Track, *, query: str | None = None, homology_protein: str | None = None,
                helixer_protein=None) -> Refinement:
    """Apply the rule to one homology gene. `gene` needs gene_name, role, contig, start, end and strand.

    `query` is the curated reference protein, `homology_protein` the translation of the homology model, and
    `helixer_protein` a function from a helixer TrackGene to its protein. With all three the extension is
    checked against the query (see the module docstring)."""
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
    start, end = min(gene.start, hx.start), max(gene.end, hx.end)
    if (start, end) == (gene.start, gene.end):
        return kept("ends_agree" if (hx.start, hx.end) == (gene.start, gene.end) else "helixer_not_longer", hx)
    verified, gain = False, None
    if query and homology_protein and helixer_protein is not None:
        verified = True
        gain = query_score(query, helixer_protein(hx)) - query_score(query, homology_protein)
        if gain < MIN_SCORE_GAIN:
            return Refinement(gene.gene_name, gene.role, gene.contig, gene.strand, KEPT,
                              "extension_not_supported_by_query", gene.start, gene.end, gene.start, gene.end,
                              hx.name, verified=True, score_gain=gain)
    d5, d3 = _deltas(gene.strand, gene.start, gene.end, start, end)
    return Refinement(gene.gene_name, gene.role, gene.contig, gene.strand, REFINED, "helixer_extension",
                      start, end, gene.start, gene.end, hx.name, d5, d3, verified, gain)


def refine_genes(genes, helixer: Track) -> list[Refinement]:
    """The rule over every present gene of a detection; one Refinement per gene, in input order."""
    return [refine_gene(g, helixer) for g in genes]
