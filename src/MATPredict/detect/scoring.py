"""Per-family fractional scoring of a gene cluster, with cross-family ambiguity detection."""
from __future__ import annotations

from dataclasses import dataclass

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey, expected_genes_for_idiomorph


@dataclass(frozen=True)
class FamilyScore:
    family_key: FamilyKey
    fraction_found: float
    genes_found: list[str]
    genes_missing: list[str]


def score_cluster(cluster: GeneCluster, families: list[Family]) -> list[FamilyScore]:
    """Score is the fraction of a family's *distinct expected gene names*
    found in the cluster -- never a summed bitscore, which would favor
    families with more genes regardless of correctness."""
    # Partition hit gene names by the family they were actually found under --
    # never pool them globally, or an identically-named gene in an unrelated
    # family (e.g. "pheromone" in both a PR family and a B-locus family)
    # would silently credit that unrelated family too.
    hit_genes_by_family: dict[FamilyKey, set[str]] = {}
    for hit in cluster.hits:
        hit_genes_by_family.setdefault(hit.family_key, set()).add(hit.gene_name)

    scores = []
    for family in families:
        found = hit_genes_by_family.get(family.key)
        if not found:
            continue
        expected = [g["name"] for g in expected_genes_for_idiomorph(family, found)]
        genes_found = [g for g in expected if g in found]
        genes_missing = [g for g in expected if g not in found]
        scores.append(FamilyScore(
            family_key=family.key,
            fraction_found=len(genes_found) / len(expected),
            genes_found=genes_found,
            genes_missing=genes_missing,
        ))
    scores.sort(key=lambda s: s.fraction_found, reverse=True)
    return scores


def is_ambiguous(scores: list[FamilyScore], floor: float = 0.5) -> bool:
    """True when 2 or more distinct families clear the floor fraction."""
    return sum(1 for s in scores if s.fraction_found >= floor) >= 2
