from __future__ import annotations
from pathlib import Path

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.search import SearchHit
from MATPredict.detect.pipeline import run_pipeline

FAMILY = Family(FamilyKey("P", "aLocus"), "pattern", None, "^a[0-9]+$",
                [{"name": "mfa1", "role": "core_MAT"}, {"name": "pra1", "role": "core_MAT"}], [1])


def test_run_pipeline_end_to_end_with_stubbed_search(tmp_path, monkeypatch):
    (tmp_path / "P").mkdir()
    (tmp_path / "P" / "order.yml").write_text(
        "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
    )

    def fake_fast_path(proteome_fasta, families, reference_fasta, runner=None):
        return [SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
                SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    def fake_genomic(*args, **kwargs):
        return []

    results = run_pipeline(
        genome_fasta=tmp_path / "genome.fa",
        proteome_fasta=tmp_path / "proteome.faa",
        taxid=None,
        db_root=tmp_path,
        reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
        search_genomic=fake_genomic,
    )
    assert len(results) == 1
    result = results[0]
    assert result.family_key == FAMILY.key
    assert result.confidence == "high"
    assert result.genes_missing == []


def test_run_pipeline_triggers_genomic_second_pass_on_missing_core_gene(tmp_path):
    (tmp_path / "P").mkdir()
    (tmp_path / "P" / "order.yml").write_text(
        "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
    )
    genomic_calls = []

    def fake_fast_path(proteome_fasta, families, reference_fasta, runner=None):
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    def fake_genomic(genome_fasta, families, reference_fasta, relaxed=False, window=None, runner=None):
        genomic_calls.append(relaxed)
        return [SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 60.0, "rec1", "exonerate_genome_relaxed")]

    results = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=fake_genomic,
    )
    assert True in genomic_calls  # relaxed second pass was actually invoked
    assert results[0].confidence == "medium"  # second-pass-confirmed core gene caps at medium
    assert results[0].genes_missing == []
