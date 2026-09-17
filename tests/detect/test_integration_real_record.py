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

    def stub_fast_path(proteome_fasta, families, reference_fasta, runner=None):
        return [
            SearchHit(f.key, gene["name"], gene["role"], "c1", 1, 100, "+", 100.0, "5270_521_aLocus_a1", "diamond_proteome")
            for f in families for gene in f.genes if f.key.locus_name == "aLocus"
        ]

    results = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=5270,
        db_root=DB_ROOT, reference_fasta=reference_fasta, search_fast_path=stub_fast_path,
        search_genomic=lambda *a, **k: [],
    )
    a_locus_result = next(r for r in results if r.family_key.locus_name == "aLocus")
    assert a_locus_result.confidence == "high"
    assert "mfa1" in a_locus_result.genes_found
    assert "pra1" in a_locus_result.genes_found
