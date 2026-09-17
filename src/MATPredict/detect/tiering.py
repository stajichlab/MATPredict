"""Per-family confidence tiering -- see spec section
"Boundary calling and confidence tiering" for the rule this encodes."""
from __future__ import annotations

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family
from MATPredict.detect.scoring import FamilyScore

_TIER_DOWNGRADE = {"high": "medium", "medium": "low", "low": "low"}


def has_flanking_conserved(family: Family) -> bool:
    return any(g["role"] == "flanking_conserved" for g in family.genes)


def assign_tier(
    score: FamilyScore,
    family: Family,
    cluster: GeneCluster,
    second_pass_used: bool,
    fragmented: bool,
) -> str:
    core_genes = {g["name"] for g in family.genes if g["role"] == "core_MAT"}
    core_found = core_genes.issubset(set(score.genes_found))

    if not core_found:
        tier = "low" if score.fraction_found == 0 else "medium"
    elif second_pass_used:
        tier = "medium"
    else:
        tier = "high"

    if fragmented and tier != "low":
        tier = _TIER_DOWNGRADE[tier]
    return tier
