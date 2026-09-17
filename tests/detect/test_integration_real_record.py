"""Integration check: the pipeline recovers the real, accepted
Mycosarcoma maydis a-locus record's genes when given its own proteins
as a stand-in 'genome proteome' and a stubbed search that returns hits
for every gene in the routed aLocus family (no live diamond/exonerate
binary required)."""
from __future__ import annotations
from pathlib import Path

import pytest

from MATPredict.detect.family_registry import load_all_families, route
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.reference_fasta import build_reference_fasta
from MATPredict.detect.search import SearchHit

DB_ROOT = Path("db")

# Confirmed real path: this record lives under "Ustilaginales" (an order),
# not "Ustilaginaceae" (a family) as the originating plan text guessed --
# see `find db/Basidiomycota -iname "5270_*"`.
REAL_RECORD = DB_ROOT / "Basidiomycota" / "Ustilaginales" / "5270_521_aLocus_a1"


@pytest.mark.skipif(not REAL_RECORD.exists(),
                     reason="requires the real curated record to be present")
def test_pipeline_recovers_real_alocus_record(tmp_path):
    families = route(5270, load_all_families(DB_ROOT))
    assert any(f.key.locus_name == "aLocus" for f in families)

    reference_fasta = build_reference_fasta(DB_ROOT, tmp_path / "reference.faa")

    def stub_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [
            SearchHit(f.key, gene["name"], gene["role"], "c1", 1, 100, "+", 100.0,
                      "5270_521_aLocus_a1", "diamond_proteome", coverage=100.0)
            for f in families for gene in f.genes if f.key.locus_name == "aLocus"
        ]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=5270,
        db_root=DB_ROOT, reference_fasta=reference_fasta, search_fast_path=stub_fast_path,
        # every gene is already found by the stubbed fast path, so no polish
        # rescue is expected -- stubbed out so no real binary can be invoked.
        polish_with_exonerate=lambda **k: None, polish_with_miniprot=lambda **k: None,
    )
    a_locus_result = next(r for r in outcome.results if r.family_key.locus_name == "aLocus")
    assert a_locus_result.confidence == "high"
    assert "mfa1" in a_locus_result.genes_found
    assert "pra1" in a_locus_result.genes_found
    # per-gene evidence survives to the result (spec section 8)
    assert {e.gene_name for e in a_locus_result.gene_evidence} >= {"mfa1", "pra1"}
    assert a_locus_result.reference_records == ["5270_521_aLocus_a1"]


@pytest.mark.skipif(not REAL_RECORD.exists(),
                     reason="requires the real curated record to be present")
def test_real_db_record_family_index_is_unambiguous():
    """Every curated record maps to exactly one (phylum, locus_name) family --
    the property that makes record-keyed hit attribution correct even where
    gene names collide across families."""
    from MATPredict.detect.family_registry import load_record_families

    record_families = load_record_families(DB_ROOT)
    families = {f.key for f in load_all_families(DB_ROOT)}
    assert record_families["5270_521_aLocus_a1"].locus_name == "aLocus"
    # every record's family is a real declared family in order.yml
    unknown = {k for k in record_families.values() if k not in families}
    assert unknown == set()
