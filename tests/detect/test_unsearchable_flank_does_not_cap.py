"""A flanking gene the run could not search must not cap the tier.

`assign_tier` already excuses `genes_not_searchable` from the CORE
requirement, and 2026-09-21 excused `optional` genes too. The flanking branch
was never given the same treatment: it asked "does this family declare a
flanking_conserved gene, and was it found?" without asking whether the run
could have found it.

Curator's ruling, 2026-09-21: "I would want to update our search model to
accept there is no flanking gene search for this clade (eg only core_MAT)."
Saccharomyces has no usable flanking GENE -- S. cerevisiae's three cassettes
(HML/MAT/HMR) are distinguished by flanking DNA, not by neighbouring genes,
and the two flanks previously declared for MATsc were measured wrong on the
S288C reference: SLA2 is not on chromosome III at all, and CHA1 abuts HML
rather than MAT. Both are now `exclude_from_search`, which puts them in
`genes_not_searchable`.

Without this fix a perfect MATsc call -- both key genes at 100% identity --
would be capped at medium for lacking a flank the run never looked for.
"""
from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.scoring import FamilyScore
from MATPredict.detect.search import SearchHit
from MATPredict.detect.tiering import assign_tier

KEY = FamilyKey("Ascomycota", "MATsc")
FAM = Family(
    key=KEY, vocabulary_type="enum", idiomorph_values=["MATa", "MATalpha"],
    idiomorph_pattern=None,
    genes=[
        {"name": "MATALPHA1", "role": "core_MAT", "present_in_idiomorphs": ["MATalpha"]},
        {"name": "MATALPHA2", "role": "core_MAT", "present_in_idiomorphs": ["MATalpha"]},
        {"name": "sla2", "role": "flanking_conserved"},
    ],
    taxonomic_scope=[4893],
)


def _cluster(*names):
    hits = [SearchHit(KEY, n, "core_MAT", "c1", 100 + i * 700, 600 + i * 700, "+",
                      100.0, "r1", "exonerate_refine") for i, n in enumerate(names)]
    return GeneCluster("c1", 100, 100 + len(names) * 700, hits)


def test_an_unsearchable_flank_does_not_cap_a_complete_core():
    found = ["MATALPHA1", "MATALPHA2"]
    score = FamilyScore(family_key=KEY, fraction_found=1.0, genes_found=found,
                        genes_missing=[], genes_not_searchable=["sla2"])
    assert assign_tier(score, FAM, _cluster(*found),
                       any_gene_unpolished=False, fragmented=False) == "high"


def test_a_SEARCHABLE_flank_that_was_not_found_still_caps():
    """The guard stays narrow: a flank the run really looked for and did not
    find is still real evidence of absence."""
    found = ["MATALPHA1", "MATALPHA2"]
    score = FamilyScore(family_key=KEY, fraction_found=0.667, genes_found=found,
                        genes_missing=["sla2"], genes_not_searchable=[])
    assert assign_tier(score, FAM, _cluster(*found),
                       any_gene_unpolished=False, fragmented=False) == "medium"


def test_a_found_flank_still_gives_high():
    found = ["MATALPHA1", "MATALPHA2", "sla2"]
    score = FamilyScore(family_key=KEY, fraction_found=1.0, genes_found=found,
                        genes_missing=[], genes_not_searchable=[])
    assert assign_tier(score, FAM, _cluster(*found),
                       any_gene_unpolished=False, fragmented=False) == "high"
