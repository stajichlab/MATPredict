# tests/detect/test_superseded_consumers.py
"""Every consumer that counts DISTINCT GENES must skip a superseded hit.

Annotating the loser of an idiomorph resolution instead of deleting it keeps
the ambiguity in the report, but it only fixes anything if the consumers
downstream stop counting it. Two failures follow from missing one:

* scoring counts one gene as two, and `expected_genes_for_idiomorph` sees two
  idiomorphs indicated so it never narrows the roster -- the inflated
  denominator this branch exists to fix.
* the evidence floor admits a cluster to polishing on `>=2 distinct genes
  including >=1 core_MAT` when there is really only one gene there.
"""
from __future__ import annotations

import dataclasses

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.idiomorph import assign_idiomorph, resolve_idiomorph_overlaps
from MATPredict.detect.pipeline import EvidenceFloor, _families_meeting_evidence_floor
from MATPredict.detect.scoring import score_cluster
from MATPredict.detect.search import SearchHit

FAM = Family(
    FamilyKey("Mucoromycota", "MAT"), "enum", ["Plus", "Minus"], None,
    [
        {"name": "tptA", "role": "flanking_conserved"},
        {"name": "sexP", "role": "core_MAT", "present_in_idiomorphs": ["Plus"]},
        {"name": "sexM", "role": "core_MAT", "present_in_idiomorphs": ["Minus"]},
        {"name": "rnhA", "role": "flanking_conserved"},
    ],
    [4827],
)

#: Every gene in FAM has a curated reference protein.
SEARCHABLE = {FAM.key: {"tptA", "sexP", "sexM", "rnhA"}}


def _hit(gene, start, end, identity, role="core_MAT"):
    return SearchHit(
        FAM.key, gene, role, "c1", start, end, "+", identity,
        f"rec_{gene}", "diamond_proteome",
    )


def _absidia_hits():
    """The real Absidia cuneospora RSA 623 locus: ground truth Plus, and the
    gene set every one of the 23 ground-truth genomes produced."""
    return [
        _hit("tptA", 14695, 15974, 75.9, role="flanking_conserved"),
        _hit("sexP", 16939, 17547, 47.3),
        _hit("sexM", 17005, 17253, 31.33),
        _hit("rnhA", 17765, 22694, 51.9, role="flanking_conserved"),
    ]


def test_scoring_ignores_a_superseded_hit_and_narrows_to_one_idiomorph():
    resolved, _ = resolve_idiomorph_overlaps(_absidia_hits(), FAM)
    cluster = GeneCluster("c1", 14695, 22694, resolved)
    score = score_cluster(cluster, [FAM], searchable_genes=SEARCHABLE)[0]
    assert score.genes_found == ["tptA", "sexP", "rnhA"]
    # Narrowed to the Plus roster, so sexM is not expected and not missing.
    assert "sexM" not in score.genes_missing
    assert score.fraction_found == 1.0


def test_the_idiomorph_is_callable_once_the_overlap_is_resolved():
    # This is the 0/23 -> 23/23 result: before resolution both idiomorphs are
    # indicated and the call is undetermined.
    hits = _absidia_hits()
    assert assign_idiomorph(FAM, [h.gene_name for h in hits]) == "undetermined"
    resolved, _ = resolve_idiomorph_overlaps(hits, FAM)
    live = [h.gene_name for h in resolved if h.superseded_by is None]
    assert assign_idiomorph(FAM, live) == "Plus"


def test_the_evidence_floor_does_not_count_a_superseded_hit_as_a_second_gene():
    # A cluster whose ONLY hits are sexM and sexP on one protein is a single
    # HMG gene. Before resolution it clears the >=2-distinct-genes bar and is
    # admitted to polishing; it must not.
    hits = [_hit("sexP", 16939, 17547, 47.3), _hit("sexM", 17005, 17253, 31.33)]
    resolved, _ = resolve_idiomorph_overlaps(hits, FAM)
    cluster = GeneCluster("c1", 16939, 17547, resolved)
    admitted = _families_meeting_evidence_floor(
        cluster, [FAM], EvidenceFloor(min_hits=2, require_core_role=True)
    )
    assert admitted == []


def test_the_evidence_floor_still_admits_a_genuine_two_gene_cluster():
    # The exclusion must bite only on the superseded hit. A real core gene
    # plus a real flank is still two distinct genes.
    resolved, _ = resolve_idiomorph_overlaps(_absidia_hits(), FAM)
    cluster = GeneCluster("c1", 14695, 22694, resolved)
    admitted = _families_meeting_evidence_floor(
        cluster, [FAM], EvidenceFloor(min_hits=2, require_core_role=True)
    )
    assert [f.key for f in admitted] == [FAM.key]


def test_a_superseded_core_hit_does_not_satisfy_the_core_role_requirement():
    # If the only core_MAT evidence in a cluster is a hit that lost an
    # idiomorph resolution, the surviving winner is still core -- but a
    # cluster where the superseded hit were the ONLY core gene must not pass
    # on it. Built directly so the winner is a flank.
    hits = [
        _hit("tptA", 100, 400, 75.9, role="flanking_conserved"),
        dataclasses.replace(_hit("sexM", 120, 380, 31.3), superseded_by="sexP"),
    ]
    cluster = GeneCluster("c1", 100, 400, hits)
    admitted = _families_meeting_evidence_floor(
        cluster, [FAM], EvidenceFloor(min_hits=1, require_core_role=True)
    )
    assert admitted == []
