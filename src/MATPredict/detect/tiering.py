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
    # A gene the run held no reference protein for is dropped from the core
    # requirement, exactly as `scoring.score_cluster` drops it from the
    # `fraction_found` denominator. Counting it here would state that the
    # genome failed to show a gene nothing could ever have found -- evidence
    # of absence manufactured from a gap in the reference database.
    #
    # Measured on Schizophyllum commune H4-8 (GCF_000143185.2, genome-only,
    # 2026-09-20): Bbeta was recovered completely -- all 8 curated genes,
    # `genes_missing=[]`, `fraction_found=1.0` -- and was still capped at
    # Medium, because the Bbeta roster's `pheromone_receptor` alias has no
    # reference protein anywhere in `db/`. Balpha (3/3 found) was demoted the
    # same way. Those were the two most complete calls in that run.
    expected_core = {
        g["name"]
        for g in expected_genes_for_idiomorph(family, score.genes_found)
        if g["role"] == "core_MAT"
    }
    core_genes = expected_core - set(score.genes_not_searchable)
    core_requirement_relaxed = core_genes != expected_core
    core_found = core_genes.issubset(set(score.genes_found))

    # Spec: Low is "a single gene hit with no other expected genes from the
    # same family found nearby". Keying Low solely on fraction_found == 0
    # made the tier unreachable in practice, because score_cluster never
    # emits a FamilyScore for a family with zero hits -- so a lone isolated
    # hit was indistinguishable from a substantial partial match. An
    # isolated single hit in a family that expects more than one gene is
    # therefore Low; any richer partial match stays Medium.
    isolated_single_hit = len(score.genes_found) <= 1 and len(family.genes) > 1

    if not core_found:
        tier = "low" if score.fraction_found == 0 or isolated_single_hit else "medium"
    elif isolated_single_hit and core_requirement_relaxed:
        # `core_found` can now be satisfied by ONE gene, when every other core
        # gene of the family is unsearchable. Without this branch that single
        # hit would go straight to High in any family with no
        # `flanking_conserved` gene -- which is every Basidiomycota family.
        # One gene is not a locus, however complete the searchable roster
        # technically was.
        #
        # Guarded on `core_requirement_relaxed` so this only ever fires where
        # the relaxation above put it: a lone core gene in a family whose
        # roster was NOT relaxed keeps its existing tier, which for a flanked
        # family with its flank missing is Medium.
        tier = "low"
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
