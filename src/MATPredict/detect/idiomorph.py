"""Idiomorph assignment, and resolution of genes that two idiomorphs both hit."""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family
from MATPredict.detect.search import SearchHit

DEFAULT_MIN_OVERLAP_FRACTION = 0.8
"""How much two hits must overlap before they are judged to be one gene.

Measured as the shared span over the SHORTER hit's span.

CALIBRATED 2026-09-20, raised from a provisional 0.5, against the corpus the
23-genome ground-truth re-run produced (236 resolution events over 23 genomes,
23 of those events on the reported loci themselves):

* across the 23 REAL loci the minimum overlap is **0.9242**, median 1.0, so
  0.8 costs nothing real and leaves 0.12 of headroom for a ragged HSP in a
  divergent taxon -- the failure mode that matters, since every one of the 23
  is Mucoromycota;
* 0.9 was equally free on this data (still 0 of 23 lost) and deliberately NOT
  taken: it leaves 0.024 of headroom on a sample of 23 genomes from one
  phylum, which is calibrating to the edge of the sample;
* raising to 0.8 also stops 4.2% of all events resolving. Those are on
  spurious clusters, and suppressing them is a mild LOSS, not a gain: an
  unresolved spurious pair keeps counting one gene as two and so inflates
  that cluster's evidence-floor gene count. This is the real argument against
  going higher still.

The two error directions remain far apart either way. Two genuinely distinct
neighbouring genes overlap at or near 0%, and where they do share bases --
normal at MAT loci, where intergenic erosion is expected -- it is a short
stretch, not four fifths of the shorter gene.

Revise from the `idiomorph_resolution` rows in the evidence diagnostics, which
carry every event's overlap fraction, not by intuition.
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


#: Methods that mean "this gene is in the annotated proteome". `search_fast_path`
#: matches the proteome with diamond; the polish stage then rewrites a hit's
#: method to the refiner's, so a gene first found by diamond and later refined
#: must not lose its proteome standing.
_PROTEOME_METHODS = frozenset(
    {"diamond_proteome", "exonerate_refine", "miniprot_refine"}
)


def _rank(hit: SearchHit) -> tuple[int, float]:
    """How good a claim this hit has to be the real gene: proteome first, then identity.

    WHICH SEARCH FOUND IT outranks identity, because the two searches mean
    different things. The fast path matches the annotated proteome, so a gene
    the annotation actually predicted is found there. The localization rescue
    then runs tblastn genome-wide for core genes MISSING from a cluster, which
    lands the other idiomorph's reference on that same locus gene precisely
    BECAUSE that idiomorph's gene is not really present. The proteome hit is
    the gene; the tblastn hit is the cross-match.

    Identity alone cannot see this, and gets it wrong as soon as a reference
    closer to the WRONG idiomorph is added: measured on the 23 ground-truth
    genomes, adding three published references pushed four Cunninghamella
    Minus genomes to Plus, dropping identity from 23/23 to 19/23, while this
    rule holds 22/23 on the same data and 16/16 on the smaller reference set.

    When both hits come from the same kind of search the path says nothing and
    identity decides, which is the honest fallback -- and the single genome
    this rule still misses (Cunninghamella polymorpha NRRL 1395) is exactly
    that case, both rescued by tblastn because the annotation predicted
    neither gene. That is the known annotation-gap problem for small MAT
    genes, not a flaw in the ranking.
    """
    return (1 if hit.method in _PROTEOME_METHODS else 0, hit.identity)


LOCUS_CLASS_MAT = "mat_locus"
LOCUS_CLASS_HOMOTHALLIC = "homothallic_candidate"
LOCUS_CLASS_IDIOMORPH_ONLY = "idiomorph_gene_only"


def classify_locus(cluster: GeneCluster, family: Family) -> str:
    """What KIND of thing this cluster is, independent of how it was admitted.

    Curator's ruling, 2026-09-20: a sub-threshold call is not junk to be thrown
    away, it is something to be filed correctly.

    * `homothallic_candidate` -- both idiomorphs' core genes present, on one
      contig, within `family.max_homothallic_separation_bp`. This is the real
      architecture of a homothallic Mucorale: Syzygites megalocarpus encodes
      both HMG transcription factors, each with its own flank. Flanking genes
      alongside make the call stronger, never weaker, so they do not change
      the class.
    * `idiomorph_gene_only` -- idiomorph-restricted core genes and NOTHING
      else: no flanking gene at all. Deliberately kept rather than discarded:
      a lone sexM or sexP is training material for a per-idiomorph HMM, which
      is a search strategy this project intends to build. It is simply not a
      locus call, so it gets its own category.
    * `mat_locus` -- everything else, i.e. a core gene with at least one
      flanking gene, the ordinary heterothallic case.

    Superseded hits are ignored throughout. A resolved cross-hit is ONE gene
    seen twice; counting it would label every ordinary heterothallic locus
    whose sexM/sexP overlap was collapsed as homothallic.
    """
    live = [
        h for h in cluster.hits
        if h.family_key == family.key and h.superseded_by is None
    ]
    if not live:
        return LOCUS_CLASS_MAT

    idiomorph_of = {g["name"]: frozenset(g.get("present_in_idiomorphs") or ())
                    for g in family.genes}

    # ONE representative hit per gene, the same one the report shows as that
    # gene's evidence. A cluster routinely holds 15+ tblastn HSPs of the same
    # region, so comparing every pairwise combination means SOME sexP/sexM pair
    # falls inside any threshold and the class stops discriminating: the
    # 44-genus sweep produced 77 homothallic candidates across 29 genera, with
    # separations up to 43,879 bp under a 20 kb bar, because a stray HSP
    # happened to sit near the other gene. Worse, the class was then decided by
    # hits the report never displayed, so the label and the coordinates
    # disagreed.
    best_by_gene: dict[str, SearchHit] = {}
    for hit in live:
        current = best_by_gene.get(hit.gene_name)
        if current is None or _rank(hit) > _rank(current):
            best_by_gene[hit.gene_name] = hit

    restricted = [h for n, h in best_by_gene.items() if idiomorph_of.get(n)]
    unrestricted = [h for n, h in best_by_gene.items() if not idiomorph_of.get(n)]

    # Both idiomorphs present and close enough to be one locus?
    for i, a in enumerate(restricted):
        for b in restricted[i + 1:]:
            if idiomorph_of[a.gene_name] & idiomorph_of[b.gene_name]:
                continue  # same idiomorph; says nothing about homothallism
            if a.contig != b.contig:
                continue
            separation = max(a.start, b.start) - min(a.end, b.end)
            if separation <= family.max_homothallic_separation_bp:
                return LOCUS_CLASS_HOMOTHALLIC

    if restricted and not unrestricted:
        return LOCUS_CLASS_IDIOMORPH_ONLY
    return LOCUS_CLASS_MAT


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
            best_a = max(group_a, key=_rank)
            best_b = max(group_b, key=_rank)
            if _rank(best_a) == _rank(best_b):
                continue  # nothing to choose between them; leave undetermined
            winner, loser = (
                (best_a, best_b) if _rank(best_a) > _rank(best_b) else (best_b, best_a)
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
