"""Per-family confidence tiering -- see spec section
"Boundary calling and confidence tiering" for the rule this encodes. When
any gene in a family is left unpolished (neither miniprot nor exonerate
--refine could confirm it), the tier is capped at Medium regardless of
flanking genes."""
from __future__ import annotations

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, expected_genes_for_idiomorph
from MATPredict.detect.scoring import FamilyScore

_TIER_DOWNGRADE = {"high": "medium", "medium": "low", "low": "low"}


def has_flanking_conserved(family: Family) -> bool:
    return any(g["role"] == "flanking_conserved" for g in family.genes)


def assign_tier(
    score: FamilyScore,
    family: Family,
    cluster: GeneCluster,
    any_gene_unpolished: bool,
    fragmented: bool,
) -> str:
    core_genes = {
        g["name"]
        for g in expected_genes_for_idiomorph(family, score.genes_found)
        if g["role"] == "core_MAT"
    }
    core_found = core_genes.issubset(set(score.genes_found))

    if not core_found:
        # Spec: Low is "a single gene hit with no other expected genes from the
        # same family found nearby". Keying Low solely on fraction_found == 0
        # made the tier unreachable in practice, because score_cluster never
        # emits a FamilyScore for a family with zero hits -- so a lone isolated
        # hit was indistinguishable from a substantial partial match. An
        # isolated single hit in a family that expects more than one gene is
        # therefore Low; any richer partial match stays Medium.
        isolated_single_hit = len(score.genes_found) <= 1 and len(family.genes) > 1
        tier = "low" if score.fraction_found == 0 or isolated_single_hit else "medium"
    elif any_gene_unpolished:
        tier = "medium"
    elif has_flanking_conserved(family):
        flanking_found = any(
            g["role"] == "flanking_conserved" and g["name"] in score.genes_found
            for g in family.genes
        )
        tier = "high" if flanking_found else "medium"
    else:
        tier = "high"

    if fragmented and tier != "low":
        tier = _TIER_DOWNGRADE[tier]
    return tier
