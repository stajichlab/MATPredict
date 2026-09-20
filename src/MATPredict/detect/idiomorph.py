"""Idiomorph assignment, and resolution of genes that two idiomorphs both hit."""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass

from MATPredict.detect.family_registry import Family
from MATPredict.detect.search import SearchHit

DEFAULT_MIN_OVERLAP_FRACTION = 0.5
"""How much two hits must overlap before they are judged to be one gene.

Measured as the shared span over the SHORTER hit's span. PROVISIONAL: the
curator ruled on 2026-09-20 that 0.5 stands "for now", to be raised once
enough resolution events accumulate in the evidence diagnostics to say what
the value should be.

0.5 was chosen because the two error directions are far apart. Across the 23
ground-truth Mucoromycota genomes every real sexM/sexP pair overlapped at
92-100%, so the bar sits well below anything observed and a ragged HSP in a
divergent taxon still resolves. Two genuinely distinct neighbouring genes, by
contrast, overlap at or near 0% -- and where they do share bases, which is
normal at MAT loci, it is a short stretch of an eroded intergenic region, not
half of the shorter gene.
"""


@dataclass(frozen=True)
class IdiomorphResolution:
    """One overlapping mutually-exclusive pair collapsed to a single gene.

    Carries both members' measurements, not just the verdict, because this
    record IS the calibration dataset for `DEFAULT_MIN_OVERLAP_FRACTION` and
    for the open question of whether coverage discriminates better than
    identity. The artifact being resolved is a shared protein DOMAIN, so a
    cross-hit should cover only part of its reference while the true gene
    covers all of it -- which makes coverage the biologically motivated
    metric and identity a proxy that happened to score 23/23. Coverage is
    recorded but not yet used; it is `None` on the tblastn and exonerate
    paths, which is itself a reason it cannot be the rule today.
    """

    contig: str
    winner: str
    loser: str
    winner_identity: float
    loser_identity: float
    overlap_fraction: float
    winner_coverage: float | None = None
    loser_coverage: float | None = None

    @property
    def margin(self) -> float:
        """Identity points separating the two. Thin margins are real: across
        the ground-truth set Plus calls are separated by 44-64 points but
        Minus calls by only 2.3-13.0, and the narrowest call seen anywhere
        (Blakeslea trispora, outside the truth set) is 0.59."""
        return self.winner_identity - self.loser_identity


def assign_idiomorph(family: Family, genes_found: list[str]) -> str:
    if family.vocabulary_type != "enum":
        return "undetermined"  # pattern (multiallelic) families: allele number isn't callable by homology alone

    indicated: set[str] = set()
    for gene in family.genes:
        if gene["name"] in genes_found:
            indicated.update(gene.get("present_in_idiomorphs", []))

    if len(indicated) == 1:
        return next(iter(indicated))
    return "undetermined"


def _idiomorphs_of(family: Family, gene_name: str) -> frozenset[str]:
    for gene in family.genes:
        if gene["name"] == gene_name:
            return frozenset(gene.get("present_in_idiomorphs") or ())
    return frozenset()


def _mutually_exclusive(family: Family, gene_a: str, gene_b: str) -> bool:
    """True when the two genes cannot belong to the same idiomorph.

    Both must name at least one idiomorph -- a gene with no
    `present_in_idiomorphs` (a flanking gene) applies to every idiomorph and
    so excludes nothing -- and their idiomorph sets must not intersect.
    """
    if gene_a == gene_b:
        return False  # two copies of one gene: duplication, normal at MAT loci
    a, b = _idiomorphs_of(family, gene_a), _idiomorphs_of(family, gene_b)
    return bool(a) and bool(b) and not (a & b)


def _overlap_fraction(a: SearchHit, b: SearchHit) -> float:
    """Shared span over the SHORTER hit's span, 0.0 on different contigs.

    Against the shorter hit, not the longer and not the union, because the
    cross-hit is typically the shorter of the two: it matches only the shared
    domain. In Absidia cuneospora the sexM cross-hit spans 249 bp inside a
    609 bp sexP hit -- 100% of the shorter but only 41% of the longer, so a
    longer-hit or reciprocal-both measure would leave the real case unresolved.
    """
    if a.contig != b.contig:
        return 0.0
    shared = min(a.end, b.end) - max(a.start, b.start) + 1
    if shared <= 0:
        return 0.0
    shorter = min(a.end - a.start + 1, b.end - b.start + 1)
    return shared / shorter


def _mutually_exclusive_pairs(
    family: Family, hits: list[SearchHit]
) -> list[tuple[str, str]]:
    """Unordered pairs of gene names present in `hits` that exclude each other."""
    names = sorted({h.gene_name for h in hits})
    return [
        (a, b)
        for i, a in enumerate(names)
        for b in names[i + 1:]
        if _mutually_exclusive(family, a, b)
    ]


def _overlap_groups(
    a_hits: list[SearchHit], b_hits: list[SearchHit], min_overlap_fraction: float
) -> list[tuple[list[SearchHit], list[SearchHit]]]:
    """Partition two gene names' hits into groups that describe one real gene.

    Two hits belong together when they overlap by at least
    `min_overlap_fraction`; grouping is transitive, so a chain of overlapping
    hits forms one group. Groups matter because gene duplication and
    multi-allele co-occurrence are normal at MAT loci: two independent copies
    of the same locus must be resolved separately, not pooled into a single
    comparison that would let a strong hit at one locus decide the call at
    the other. A group with no member from BOTH names has nothing to resolve
    and is dropped.
    """
    #: Union-find over the combined hit list, keyed by position.
    combined = a_hits + b_hits
    parent = list(range(len(combined)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(combined)):
        for j in range(i + 1, len(combined)):
            if _overlap_fraction(combined[i], combined[j]) >= min_overlap_fraction:
                parent[find(i)] = find(j)

    groups: dict[int, tuple[list[SearchHit], list[SearchHit]]] = {}
    for index, hit in enumerate(combined):
        group = groups.setdefault(find(index), ([], []))
        group[0 if index < len(a_hits) else 1].append(hit)
    return [(ga, gb) for ga, gb in groups.values() if ga and gb]


def resolve_idiomorph_overlaps(
    hits: list[SearchHit],
    family: Family,
    min_overlap_fraction: float = DEFAULT_MIN_OVERLAP_FRACTION,
) -> tuple[list[SearchHit], list[IdiomorphResolution]]:
    """Collapse pairs of mutually exclusive idiomorph genes that hit one gene.

    Returns the hits with each loser annotated `superseded_by`, plus one
    `IdiomorphResolution` per collapsed pair. Nothing is removed.

    Run this EARLY, right after clustering, so that every consumer counting
    distinct genes sees the resolved picture. Deferring it to scoring would
    fix the idiomorph call but leave the polish-admission bar inflated: a
    cluster whose only hits are sexM and sexP on one protein clears
    `EvidenceFloor(min_hits=2, require_core_role=True)` on the strength of a
    single gene.

    An exact identity tie is left unresolved. There is nothing to choose
    between the two, the choice would flip on floating-point noise, and
    reporting the locus as `undetermined` is the honest outcome.
    """
    superseded: dict[int, str] = {}
    events: list[IdiomorphResolution] = []
    for gene_a, gene_b in _mutually_exclusive_pairs(family, hits):
        a_hits = [h for h in hits if h.gene_name == gene_a]
        b_hits = [h for h in hits if h.gene_name == gene_b]
        for group_a, group_b in _overlap_groups(a_hits, b_hits, min_overlap_fraction):
            # Best against best, ONE verdict per real gene. A curated family
            # routinely contributes several reference proteins per gene name
            # (the Mucoromycota database holds 3 sexP and 3 sexM), so one
            # locus gene draws a fistful of hits under each name. Comparing
            # every hit of one name against every hit of the other yields
            # contradictory verdicts -- on the real Absidia cuneospora locus,
            # 9 of them, one with sexM beating sexP -- and a margin taken from
            # the closest accidental pairing rather than from the actual call.
            best_a = max(group_a, key=lambda h: h.identity)
            best_b = max(group_b, key=lambda h: h.identity)
            if best_a.identity == best_b.identity:
                continue  # nothing to choose between them; leave undetermined
            winner, loser = (
                (best_a, best_b) if best_a.identity > best_b.identity else (best_b, best_a)
            )
            losing_group = group_b if winner is best_a else group_a
            for hit in losing_group:
                superseded[id(hit)] = winner.gene_name
            events.append(IdiomorphResolution(
                contig=winner.contig,
                winner=winner.gene_name,
                loser=loser.gene_name,
                winner_identity=winner.identity,
                loser_identity=loser.identity,
                overlap_fraction=_overlap_fraction(winner, loser),
                winner_coverage=winner.coverage,
                loser_coverage=loser.coverage,
            ))
    resolved = [
        dataclasses.replace(h, superseded_by=superseded[id(h)]) if id(h) in superseded else h
        for h in hits
    ]
    return resolved, events
