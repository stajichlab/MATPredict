"""An `optional` flanking gene is a BONUS, never a penalty.

Curator's ruling, 2026-09-21: "refine the sla2 flanking for lineages outside of
Saccharomyces ... perhaps sla2 is a bonus flank to search for but not to
penalize if not present?"

The measured problem `sla2` poses for the `MATsc` family, which spans the whole
of Saccharomycetaceae:

* in *Saccharomyces*, SLA2 is not on chromosome III AT ALL (checked against the
  S288C reference, 150 annotated chrIII genes), so it is not a MAT flank there;
* in *Kluyveromyces lactis* and *Lachancea thermotolerans* it IS genuinely
  adjacent to MAT (sla2 -> MATA1 -> MATA2 on NC_006039.1 and NC_013082.1).

One flag cannot say "wrong here, right there". `exclude_from_search: true`
silenced it for the whole family and threw away the genera where it works.
`optional: true` is the shape that fits: searched everywhere, reported when
found, promoting a call to High when present, and costing nothing where the
architecture simply does not have it.

`optional` already keeps a gene out of `fraction_found` (scoring.py) and out of
the CORE requirement (assign_tier). This adds the missing third place: the
FLANKING requirement.
"""
from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.scoring import FamilyScore
from MATPredict.detect.search import SearchHit
from MATPredict.detect.tiering import assign_tier

KEY = FamilyKey("Ascomycota", "MATsc")


def _fam(sla2_optional: bool):
    sla2 = {"name": "sla2", "role": "flanking_conserved"}
    if sla2_optional:
        sla2["optional"] = True
    return Family(
        key=KEY, vocabulary_type="enum", idiomorph_values=["MATa", "MATalpha"],
        idiomorph_pattern=None,
        genes=[
            {"name": "MATALPHA1", "role": "core_MAT", "present_in_idiomorphs": ["MATalpha"]},
            {"name": "MATALPHA2", "role": "core_MAT", "present_in_idiomorphs": ["MATalpha"]},
            sla2,
        ],
        taxonomic_scope=[4893],
    )


def _cluster(*names):
    hits = [SearchHit(KEY, n, "core_MAT", "c1", 100 + i * 800, 700 + i * 800, "+",
                      100.0, "r1", "exonerate_refine") for i, n in enumerate(names)]
    return GeneCluster("c1", 100, 100 + len(names) * 800, hits)


def _tier(fam, found, missing):
    score = FamilyScore(family_key=KEY, fraction_found=1.0, genes_found=list(found),
                        genes_missing=list(missing), genes_not_searchable=[])
    return assign_tier(score, fam, _cluster(*found),
                       any_gene_unpolished=False, fragmented=False)


def test_a_missing_optional_flank_does_not_cap_the_tier():
    """Saccharomyces: the core is complete, sla2 was searched and is genuinely
    absent, and that must not demote the call."""
    assert _tier(_fam(True), ["MATALPHA1", "MATALPHA2"], ["sla2"]) == "high"


def test_a_found_optional_flank_still_earns_high():
    """Kluyveromyces/Lachancea: sla2 IS adjacent there, and finding it is
    corroboration, not a new requirement."""
    assert _tier(_fam(True), ["MATALPHA1", "MATALPHA2", "sla2"], []) == "high"


def test_a_REQUIRED_flank_that_is_missing_still_caps():
    """The guard stays narrow. Only `optional` is excused -- a required flank
    the run searched for and did not find is still evidence of absence."""
    assert _tier(_fam(False), ["MATALPHA1", "MATALPHA2"], ["sla2"]) == "medium"


def test_a_required_flank_that_is_found_still_earns_high():
    assert _tier(_fam(False), ["MATALPHA1", "MATALPHA2", "sla2"], []) == "high"
