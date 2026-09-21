# tests/detect/test_tiering.py
from __future__ import annotations

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.scoring import FamilyScore
from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.tiering import assign_tier, has_flanking_conserved

FLANKLESS = Family(FamilyKey("P", "aLocus"), "pattern", None, "^a[0-9]+$",
                    [{"name": "mfa1", "role": "core_MAT"}, {"name": "pra1", "role": "core_MAT"}], [1])
FLANKED = Family(FamilyKey("P", "MAT"), "enum", ["a", "alpha"], None,
                  [{"name": "STE3", "role": "core_MAT"}, {"name": "flank1", "role": "flanking_conserved"}], [1])


def test_has_flanking_conserved():
    assert has_flanking_conserved(FLANKED) is True
    assert has_flanking_conserved(FLANKLESS) is False


def test_flankless_family_high_tier_on_all_core_genes_found():
    score = FamilyScore(FLANKLESS.key, 1.0, ["mfa1", "pra1"], [])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, FLANKLESS, cluster, any_gene_unpolished=False, fragmented=False) == "high"


def test_flanked_family_medium_tier_when_core_gene_only_found_via_second_pass():
    score = FamilyScore(FLANKED.key, 1.0, ["STE3", "flank1"], [])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, FLANKED, cluster, any_gene_unpolished=True, fragmented=False) == "medium"


def test_flanked_family_medium_tier_when_no_flanking_gene_found():
    score = FamilyScore(FLANKED.key, 0.5, ["STE3"], ["flank1"])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, FLANKED, cluster, any_gene_unpolished=False, fragmented=False) == "medium"


def test_isolated_single_hit_is_low():
    """Spec's Low tier: "a single gene hit with no other expected genes from
    the same family found nearby". Previously this returned "medium" and Low
    was reachable only at fraction_found == 0, which score_cluster never
    emits -- making the tier dead code."""
    score = FamilyScore(FLANKLESS.key, 0.5, ["mfa1"], ["pra1"])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, FLANKLESS, cluster, any_gene_unpolished=False, fragmented=False) == "low"


def test_partial_match_with_several_genes_is_medium():
    """More than one gene found, but not the full core set -> Medium, not Low."""
    family = Family(FamilyKey("P", "big"), "pattern", None, "^a[0-9]+$",
                    [{"name": "g1", "role": "core_MAT"}, {"name": "g2", "role": "core_MAT"},
                     {"name": "g3", "role": "core_MAT"}, {"name": "g4", "role": "core_MAT"}], [1])
    score = FamilyScore(family.key, 0.5, ["g1", "g2"], ["g3", "g4"])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, family, cluster, any_gene_unpolished=False, fragmented=False) == "medium"


def test_fragmented_locus_downgraded_one_tier():
    score = FamilyScore(FLANKLESS.key, 1.0, ["mfa1", "pra1"], [])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, FLANKLESS, cluster, any_gene_unpolished=False, fragmented=True) == "medium"


# Family split across 2 idiomorphs, mirroring this session's real regression:
# a single-idiomorph genome only ever has half the family's core genes, so
# core_found against the FULL roster (both idiomorphs) could never be
# satisfied for a real single-idiomorph genome.
FAM_IDIOMORPHIC = Family(
    FamilyKey("P", "MAT"), "enum", ["MAT1-1", "MAT1-2"], None,
    [
        {"name": "a1", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-1"]},
        {"name": "a2", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-1"]},
        {"name": "b1", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-2"]},
        {"name": "b2", "role": "core_MAT", "present_in_idiomorphs": ["MAT1-2"]},
    ],
    [1],
)


def test_single_idiomorph_genome_reaches_high_tier():
    # Only MAT1-2's genes (b1, b2) are found -- against the full 4-gene
    # roster core_found would require a1/a2 too and could never be
    # satisfied. Narrowed to MAT1-2 alone, both core genes are found.
    score = FamilyScore(FAM_IDIOMORPHIC.key, 0.5, ["b1", "b2"], ["a1", "a2"])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, FAM_IDIOMORPHIC, cluster, any_gene_unpolished=False, fragmented=False) == "high"


def test_both_idiomorphs_found_still_requires_full_core_set():
    # Only one idiomorph's genes found (b1, b2) plus one gene (a1) from the
    # other -- found idiomorphs are {MAT1-1, MAT1-2}, so no narrowing
    # happens and core_found still requires all 4 genes (unchanged, full-
    # roster behavior for the both-idiomorphs-present case).
    score = FamilyScore(FAM_IDIOMORPHIC.key, 0.75, ["a1", "b1", "b2"], ["a2"])
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(score, FAM_IDIOMORPHIC, cluster, any_gene_unpolished=False, fragmented=False) == "medium"


# --- Unsearchable core genes must not block `core_found` -------------------
# `scoring.score_cluster` already drops genes with no reference protein from
# the `fraction_found` denominator and reports them as `genes_not_searchable`.
# `assign_tier` never read that field, so a gene nothing could ever find was
# still counted as a core requirement the genome had failed to meet.
#
# Measured on Schizophyllum commune H4-8 (GCF_000143185.2), genome-only,
# 2026-09-20: the Bbeta locus was recovered COMPLETELY -- all 8 curated genes
# found, `genes_missing=[]`, `fraction_found=1.0` -- and was still capped at
# `medium`, because the roster's `pheromone_receptor` alias has no reference
# protein anywhere in db/ and stayed in `core_genes`. Balpha (3/3 found) was
# demoted the same way. Those were the two most complete calls in the run.

BBETA_LIKE = Family(
    FamilyKey("Basidiomycota", "Bbeta"), "pattern", None, "^[0-9]+$",
    [{"name": "pheromone_receptor", "role": "core_MAT"},  # no reference protein
     {"name": "bbr2", "role": "core_MAT"},
     {"name": "bbp2-1", "role": "core_MAT"},
     {"name": "bbp2-2", "role": "core_MAT"}], [5334])


def test_unsearchable_core_gene_does_not_block_high_tier():
    """Every SEARCHABLE core gene found, nothing missing -> High.

    The one unmet core gene has no reference protein, so its absence is not
    evidence of anything about this genome."""
    score = FamilyScore(
        BBETA_LIKE.key, 1.0,
        genes_found=["bbr2", "bbp2-1", "bbp2-2"],
        genes_missing=[],
        genes_not_searchable=["pheromone_receptor"],
    )
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(
        score, BBETA_LIKE, cluster, any_gene_unpolished=False, fragmented=False
    ) == "high"


def test_unsearchable_genes_do_not_promote_a_lone_hit_to_high():
    """The complement of the above, and the reason the fix is not merely
    `core_genes -= not_searchable`: if every other core gene is unsearchable,
    ONE found gene would otherwise satisfy `core_found` and -- because no
    Basidiomycota family declares a `flanking_conserved` gene -- go straight
    to High. An isolated single hit stays Low."""
    score = FamilyScore(
        BBETA_LIKE.key, 1.0,
        genes_found=["bbr2"],
        genes_missing=[],
        genes_not_searchable=["pheromone_receptor", "bbp2-1", "bbp2-2"],
    )
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(
        score, BBETA_LIKE, cluster, any_gene_unpolished=False, fragmented=False
    ) == "low"


def test_searchable_core_gene_still_missing_is_not_high():
    """Guard against over-relaxing: a genuinely missing SEARCHABLE core gene
    must still keep the call below High."""
    score = FamilyScore(
        BBETA_LIKE.key, 0.67,
        genes_found=["bbr2", "bbp2-1"],
        genes_missing=["bbp2-2"],
        genes_not_searchable=["pheromone_receptor"],
    )
    cluster = GeneCluster("c1", 1, 100, [])
    assert assign_tier(
        score, BBETA_LIKE, cluster, any_gene_unpolished=False, fragmented=False
    ) == "medium"
