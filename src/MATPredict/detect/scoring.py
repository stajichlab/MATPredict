"""Per-family fractional scoring of a gene cluster, with cross-family ambiguity detection."""
from __future__ import annotations

from dataclasses import dataclass, field

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey, expected_genes_for_idiomorph


@dataclass(frozen=True)
class FamilyScore:
    family_key: FamilyKey
    fraction_found: float
    genes_found: list[str]
    genes_missing: list[str]
    genes_not_searchable: list[str] = field(default_factory=list)
    """Expected genes this run had no way to find, excluded from the denominator.

    A gene is here when the reference set the run actually searched with held
    no protein for it -- because no curated record carries that gene, or
    because phylum routing left its records out. Such a gene cannot be found
    by any genome, so counting it as missing would not be conservative, it
    would be meaningless: the ceiling on `fraction_found` would then track
    gaps in the curated database rather than the biology being scored.
    """
    genes_optional_found: list[str] = field(default_factory=list)
    """Optional genes that WERE found. Informative, never scored.

    An `optional: true` gene is outside `fraction_found` entirely -- out of
    the denominator AND out of the numerator -- because the biology does not
    require it. It is still searched and still reported, in `genes_found` and
    named again here, because its presence carries information even though
    its absence carries none. Curator ruling, J. Stajich, 2026-09-20:
    "btbA shouldn't count in the denominator but we should know when it is
    present."

    This is why `fraction_found` is NOT derivable from `genes_found` and
    `genes_missing` alone: subtract this list from `genes_found` first.
    """


def score_cluster(
    cluster: GeneCluster,
    families: list[Family],
    searchable_genes: dict[FamilyKey, set[str]] | None = None,
) -> list[FamilyScore]:
    """Score is the fraction of a family's *distinct expected gene names*
    found in the cluster -- never a summed bitscore, which would favor
    families with more genes regardless of correctness.

    `searchable_genes` names, per family, the genes the run's reference set
    actually held a protein for; expected genes outside it are dropped from
    the denominator and reported as `genes_not_searchable` instead. Passing
    nothing -- or a map with no entry for a family -- means "no information
    about searchability" and keeps the whole roster in the denominator, which
    is the pre-existing behaviour every current caller relies on. It never
    means "nothing is searchable", which would divide by zero.

    Searchability is applied AFTER the idiomorph roster is narrowed, so the
    denominator is the intersection of "this idiomorph expects it" and "we
    could have found it". Applied the other way round, the other idiomorph's
    genes would be reported as unsearchable.
    """
    # Partition hit gene names by the family they were actually found under --
    # never pool them globally, or an identically-named gene in an unrelated
    # family (e.g. "pheromone" in both a PR family and a B-locus family)
    # would silently credit that unrelated family too.
    hit_genes_by_family: dict[FamilyKey, set[str]] = {}
    for hit in cluster.hits:
        # A superseded hit lost an idiomorph resolution: it and the winner hit
        # the SAME locus gene, so counting it would score one gene as two and
        # would leave two idiomorphs indicated, which stops
        # `expected_genes_for_idiomorph` narrowing the roster at all.
        if hit.superseded_by is not None:
            continue
        hit_genes_by_family.setdefault(hit.family_key, set()).add(hit.gene_name)

    scores = []
    for family in families:
        found = hit_genes_by_family.get(family.key)
        if not found:
            continue
        expected_genes = expected_genes_for_idiomorph(family, found)
        # An `optional: true` gene leaves the score entirely -- numerator and
        # denominator both -- and is reported separately instead. Measured
        # reason (283-genome BFD sweep): btbA is Plus-restricted, so it sat in
        # every Plus locus's expected roster, yet only 26% of Plus loci have
        # it. Plus was scored out of 6 and Minus out of 5, so a complete Plus
        # locus with no btbA scored 5/6 = 0.833 and could never reach 1.0
        # while the equivalent Minus locus did. That is a database artefact
        # about one gene's single reference protein, not biology.
        optional = {g["name"] for g in expected_genes if g.get("optional")}
        expected = [g["name"] for g in expected_genes if g["name"] not in optional]
        family_searchable = (
            None if searchable_genes is None else searchable_genes.get(family.key)
        )
        if family_searchable is None:
            searchable, not_searchable = expected, []
        else:
            # A gene that was FOUND is searchable by definition, whatever the
            # map says. `searchable_genes` is an inference about what the run
            # could have found and it can be wrong in the permissive
            # direction: the short-ORF exclusion drops genes whose only
            # reference protein is too short to localize reliably, yet the
            # windowed polish rescue finds exactly such genes -- that rescue
            # is the reason this pipeline exists. Trusting the map over the
            # evidence would delete a real detection from `genes_found`.
            searchable = [g for g in expected if g in family_searchable or g in found]
            not_searchable = [
                g for g in expected if g not in family_searchable and g not in found
            ]
        scored_found = [g for g in searchable if g in found]
        genes_missing = [g for g in searchable if g not in found]
        # Optional genes are appended to `genes_found` so a reader sees what
        # was actually there, but they are NOT in `scored_found`, which is
        # the numerator. Order is preserved from the family roster so the
        # report reads in locus order.
        # A `slot` groups optional, idiomorph-alternative genes (SXI1 alpha /
        # SXI2 a) into ONE bonus unit: when any member is found the unit adds
        # 1 to numerator AND denominator; when none is, it stays out of the
        # score like any optional gene. So it can raise a score, never lower
        # one -- the curator ruled on 2026-09-21 that SXI1/SXI2 are not
        # established outside Cryptococcus, so their absence must cost nothing.
        # Measured need (results/2026-09-26_cryptococcus_zero/NOTE.md): four
        # MATalpha C. neoformans assemblies split the locus over 5-6 contigs;
        # the true piece (SXI1 + FAO1) scored 1/3 and fell below the 0.5 floor.
        #
        # A slot scores only when EXACTLY ONE member is found. The members
        # are alternatives, so two of them in one cluster is what a
        # homeodomain paralog hitting both references looks like, not
        # evidence of either idiomorph. Measured on 243 Cryptococcus genomes
        # (results/2026-09-26_sxi_slot_gateA/): 12 new calls held both, and 3
        # of them turned a single-idiomorph genome into a+alpha.
        members_found: dict[str, list[str]] = {}
        for g in expected_genes:
            if g.get("slot") and g["name"] in optional and g["name"] in found:
                members_found.setdefault(g["slot"], []).append(g["name"])
        slots_hit = {slot for slot, names in members_found.items() if len(names) == 1}
        slot_found = [n for slot in slots_hit for n in members_found[slot]]
        optional_found = [
            g["name"] for g in expected_genes
            if g["name"] in optional and g["name"] in found
            and g["name"] not in slot_found
        ]
        if not searchable:
            # Every expected gene was optional or unsearchable, and the only
            # hits are optional ones. There is no required gene to score, and
            # 0/0 is not a score -- an optional gene must never be able to
            # manufacture a fraction. Report it and score zero. A slot is
            # optional too, so it cannot manufacture one either.
            fraction = 0.0
        else:
            fraction = (len(scored_found) + len(slots_hit)) / (len(searchable) + len(slots_hit))
        scores.append(FamilyScore(
            family_key=family.key,
            fraction_found=fraction,
            genes_found=scored_found + slot_found + optional_found,
            genes_missing=genes_missing,
            genes_not_searchable=not_searchable,
            genes_optional_found=optional_found,
        ))
    scores.sort(key=lambda s: s.fraction_found, reverse=True)
    return scores


def is_ambiguous(scores: list[FamilyScore], floor: float = 0.5) -> bool:
    """True when 2 or more distinct families clear the floor fraction."""
    return sum(1 for s in scores if s.fraction_found >= floor) >= 2


def count_distinct_intervals(hits, min_overlap_fraction: float = 0.5) -> int:
    """How many separate places on the genome these hits actually occupy.

    Hits are grouped transitively when they overlap by at least
    `min_overlap_fraction` of the SHORTER hit -- the same test
    `idiomorph._overlap_groups` uses, so one cluster is never measured two
    different ways. Different contigs never group.

    **Why this exists.** Gene COUNT stopped discriminating once it emerged
    that several roster entries could name one ORF: on 334 Tremellales
    genomes, 927 medium-confidence calls were a single ~95 bp interval
    reported as three genes. Interval count is not fooled by that. Measured on
    the same panel after the reference dedup: every one of the 2,157 low calls
    occupies exactly ONE interval, while all 134 high and all 435 medium calls
    occupy two or more.

    **Safe for tandem duplication**, which is real at Basidiomycota MAT loci:
    tandem copies sit ADJACENT rather than overlapping, so they group
    separately and are counted separately. That is the property name-collapsing
    does not have, and the reason both changes were needed.

    Superseded hits are excluded. A resolved idiomorph cross-match is one gene
    seen twice under two names; the loser stays in the report as evidence and
    must not inflate the count.

    Reported as a field, not used as a gate. See the report writer.
    """
    live = [h for h in hits if getattr(h, "superseded_by", None) is None]
    if not live:
        return 0
    parent = list(range(len(live)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(live)):
        for j in range(i + 1, len(live)):
            a, b = live[i], live[j]
            if a.contig != b.contig:
                continue
            shared = min(a.end, b.end) - max(a.start, b.start) + 1
            if shared <= 0:
                continue
            shorter = min(a.end - a.start + 1, b.end - b.start + 1)
            if shared / shorter >= min_overlap_fraction:
                parent[find(i)] = find(j)
    return len({find(i) for i in range(len(live))})
