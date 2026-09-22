"""An `optional: true` gene must not be part of the core requirement.

`scoring.score_cluster` already drops optional genes from `fraction_found`
(numerator and denominator). `tiering.assign_tier` did not: it built
`expected_core` from every `core_MAT` gene in the roster and required ALL of
them to be found. So marking a gene optional removed it from the score but
still let its absence cap the tier -- the same "requiring something that is
not required" shape as `genes_not_searchable`, which this function already
fixes two lines above.

It became load-bearing with the 2026-09-21 Cryptococcus recuration: the
curator ruled that the recombination-trapped genes (STE20, RPO41, RPL39) and
the homeodomain genes (SXI1, SXI2) go in as `optional`, because they are
expected in Cryptococcus but not known to be universal across Tremellales.
Without this fix a real MAT locus showing the pheromone/receptor core would be
demoted for lacking genes the curator explicitly said are not required.
"""
from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.scoring import FamilyScore
from MATPredict.detect.search import SearchHit
from MATPredict.detect.tiering import assign_tier

KEY = FamilyKey("Basidiomycota", "MAT")

FAM = Family(
    key=KEY, vocabulary_type="enum", idiomorph_values=["a", "alpha"],
    idiomorph_pattern=None,
    genes=[
        {"name": "MFalpha", "role": "core_MAT", "present_in_idiomorphs": ["alpha"]},
        {"name": "STE3", "role": "core_MAT"},
        {"name": "SXI1", "role": "core_MAT", "present_in_idiomorphs": ["alpha"],
         "optional": True},
        {"name": "STE20", "role": "core_MAT", "optional": True},
        {"name": "FAO1", "role": "flanking_conserved"},
    ],
    taxonomic_scope=[5234],
)


def _cluster(*names):
    hits = [
        SearchHit(KEY, n, "core_MAT", "c1", 100 + i * 1000, 600 + i * 1000, "+",
                  95.0, "rec1", "diamond_proteome")
        for i, n in enumerate(names)
    ]
    return GeneCluster("c1", 100, 100 + len(names) * 1000, hits)


def test_a_missing_optional_core_gene_does_not_block_the_high_tier():
    found = ["MFalpha", "STE3", "FAO1"]
    score = FamilyScore(
        family_key=KEY, fraction_found=1.0, genes_found=found,
        genes_missing=[], genes_not_searchable=[],
    )
    tier = assign_tier(score, FAM, _cluster(*found),
                       any_gene_unpolished=False, fragmented=False)
    assert tier == "high", "SXI1/STE20 are optional; their absence must not demote"


def test_a_missing_REQUIRED_core_gene_still_blocks_it():
    """The guard must stay narrow -- only `optional` is excused."""
    found = ["MFalpha", "FAO1"]  # STE3 is required and absent
    score = FamilyScore(
        family_key=KEY, fraction_found=0.5, genes_found=found,
        genes_missing=["STE3"], genes_not_searchable=[],
    )
    tier = assign_tier(score, FAM, _cluster(*found),
                       any_gene_unpolished=False, fragmented=False)
    assert tier != "high"


def test_finding_an_optional_gene_still_helps_nothing_break():
    found = ["MFalpha", "STE3", "SXI1", "STE20", "FAO1"]
    score = FamilyScore(
        family_key=KEY, fraction_found=1.0, genes_found=found,
        genes_missing=[], genes_not_searchable=[],
    )
    assert assign_tier(score, FAM, _cluster(*found),
                       any_gene_unpolished=False, fragmented=False) == "high"
