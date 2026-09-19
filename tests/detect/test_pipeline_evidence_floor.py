from __future__ import annotations

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.pipeline import EvidenceFloor, _families_meeting_evidence_floor
from MATPredict.detect.search import SearchHit


def _family(phylum, name, genes=None):
    return Family(
        key=FamilyKey(phylum, name), vocabulary_type="enum",
        idiomorph_values=["a", "alpha"], idiomorph_pattern=None,
        genes=genes or [{"name": "G1", "role": "core_MAT"}],
        taxonomic_scope=[1],
    )


def _hit(family_key, gene_name="G1", role="core_MAT", identity=30.0):
    return SearchHit(
        family_key=family_key, gene_name=gene_name, role=role, contig="c1",
        start=1, end=100, strand="+", identity=identity,
        reference_record_id="rec1", method="tblastn_genome",
    )


def test_default_evidence_floor_admits_any_single_hit_regardless_of_identity():
    real = _family("Ascomycota", "MAT")
    spurious = _family("Basidiomycota", "bLocus")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(real.key, identity=34.9), _hit(spurious.key, identity=72.7),
    ])

    result = _families_meeting_evidence_floor(cluster, [real, spurious], EvidenceFloor())

    assert {f.key for f in result} == {real.key, spurious.key}


def test_evidence_floor_with_no_hits_admits_nothing():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[])

    assert _families_meeting_evidence_floor(cluster, [family], EvidenceFloor()) == []


def test_min_hits_floor_rejects_a_single_hit_family():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[_hit(family.key)])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(min_hits=2))

    assert result == []


def test_min_hits_floor_admits_a_family_with_enough_distinct_genes():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(family.key, gene_name="G1"), _hit(family.key, gene_name="G2"),
    ])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(min_hits=2))

    assert [f.key for f in result] == [family.key]


def test_min_identity_floor_rejects_a_family_whose_best_hit_is_below_it():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[_hit(family.key, identity=34.9)])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(min_identity=50.0))

    assert result == []


def test_require_core_role_rejects_a_family_with_only_flanking_hits():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(family.key, role="flanking_conserved", identity=90.0),
    ])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(require_core_role=True))

    assert result == []


def test_require_core_role_admits_a_family_with_at_least_one_core_hit():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(family.key, role="flanking_conserved", identity=90.0),
        _hit(family.key, gene_name="G1", role="core_MAT", identity=40.0),
    ])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(require_core_role=True))

    assert [f.key for f in result] == [family.key]


def test_evidence_diagnostics_written_when_path_given(tmp_path):
    from MATPredict.detect.pipeline import _write_evidence_diagnostics

    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(family.key, gene_name="G1", role="core_MAT", identity=34.9),
    ])
    out_path = tmp_path / "diagnostics.jsonl"

    _write_evidence_diagnostics(out_path, cluster, family, admitted=True)
    _write_evidence_diagnostics(out_path, cluster, family, admitted=False)

    import json
    lines = out_path.read_text().splitlines()
    assert len(lines) == 2
    row = json.loads(lines[0])
    assert row["family"] == "Ascomycota:MAT"
    assert row["contig"] == "c1"
    assert row["gene_count"] == 1
    assert row["best_identity"] == 34.9
    assert row["admitted"] is True


def test_evidence_diagnostics_hit_count_tracks_raw_hsps_separately_from_gene_count(tmp_path):
    from MATPredict.detect.pipeline import _write_evidence_diagnostics

    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(family.key, gene_name="G1", role="core_MAT", identity=34.9),
        _hit(family.key, gene_name="G1", role="core_MAT", identity=40.0),
    ])
    out_path = tmp_path / "diagnostics.jsonl"

    _write_evidence_diagnostics(out_path, cluster, family, admitted=True)

    import json
    row = json.loads(out_path.read_text().splitlines()[0])
    assert row["gene_count"] == 1
    assert row["hit_count"] == 2


def test_min_hits_now_counts_distinct_genes_not_raw_hsps():
    family = _family("Ascomycota", "MAT")
    # Two HSPs for the SAME gene (e.g. two curated reference records both hit
    # by tblastn) -- this is exactly the case search.py's own lack of
    # per-gene dedup produces, and it must NOT count as 2 toward min_hits.
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(family.key, gene_name="G1", identity=40.0),
        _hit(family.key, gene_name="G1", identity=42.0),
    ])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(min_hits=2))

    assert result == []  # only 1 DISTINCT gene, even though there are 2 hit objects


def test_min_hits_admits_when_distinct_gene_count_clears_the_floor():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(family.key, gene_name="G1", identity=40.0),
        _hit(family.key, gene_name="G1", identity=42.0),  # 2nd HSP, same gene
        _hit(family.key, gene_name="G2", identity=38.0),  # a genuinely different gene
    ])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(min_hits=2))

    assert [f.key for f in result] == [family.key]  # 2 distinct genes (G1, G2)
