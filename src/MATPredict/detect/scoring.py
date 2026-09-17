"""Per-family fractional scoring of a gene cluster, with cross-family ambiguity detection."""
from __future__ import annotations

from dataclasses import dataclass

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey


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
    # Collect all distinct gene names from hits, regardless of family association
    hit_gene_names: set[str] = set()
    for hit in cluster.hits:
        hit_gene_names.add(hit.gene_name)

    scores = []
    for family in families:
        expected = [g["name"] for g in family.genes]
        genes_found = [g for g in expected if g in hit_gene_names]
        if not genes_found:
            continue
        genes_missing = [g for g in expected if g not in hit_gene_names]
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
