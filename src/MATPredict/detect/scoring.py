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
        expected = [g["name"] for g in expected_genes_for_idiomorph(family, found)]
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
        genes_found = [g for g in searchable if g in found]
        genes_missing = [g for g in searchable if g not in found]
        scores.append(FamilyScore(
            family_key=family.key,
            # `searchable` cannot be empty here: this family has at least one
            # hit, and every found gene is kept searchable above regardless of
            # what the map claims.
            fraction_found=len(genes_found) / len(searchable),
            genes_found=genes_found,
            genes_missing=genes_missing,
            genes_not_searchable=not_searchable,
        ))
    scores.sort(key=lambda s: s.fraction_found, reverse=True)
    return scores


def is_ambiguous(scores: list[FamilyScore], floor: float = 0.5) -> bool:
    """True when 2 or more distinct families clear the floor fraction."""
    return sum(1 for s in scores if s.fraction_found >= floor) >= 2
