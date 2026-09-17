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


def test_second_pass_in_one_cluster_does_not_cap_a_different_cluster_of_the_same_family(tmp_path):
    """Same family, two independent spatial clusters (e.g. gene-duplication /
    multi-allele co-occurrence). Cluster on contig c1 needs the relaxed
    genomic second pass to confirm mfa1. Cluster on contig c2 gets both core
    genes from the fast path alone -- it must reach "high", not be capped at
    "medium" just because the c1 cluster (same family) needed a second pass."""
    (tmp_path / "P").mkdir()
    (tmp_path / "P" / "order.yml").write_text(
        "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
    )

    def fake_fast_path(proteome_fasta, families, reference_fasta, runner=None):
        return [
            # c1 cluster: only pra1 found by the fast path -> needs second pass.
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
            # c2 cluster: both core genes found by the fast path -> no second pass needed.
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c2", 100, 200, "+", 95.0, "rec2", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c2", 300, 400, "+", 95.0, "rec2", "diamond_proteome"),
        ]

    def fake_genomic(genome_fasta, families, reference_fasta, relaxed=False, window=None, runner=None):
        contig = window[0] if window else None
        if contig == "c1":
            return [SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 60.0, "rec1",
                               "exonerate_genome_relaxed")]
        return []

    results = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=fake_genomic,
    )

    by_contig = {r.contig: r for r in results}
    assert len(results) == 2
    assert by_contig["c1"].confidence == "medium"  # this cluster's own second pass caps it
    assert by_contig["c2"].confidence == "high"  # unrelated cluster, same family, must NOT be capped


def test_short_orf_gene_reported_as_not_searchable_not_missing(tmp_path):
    (tmp_path / "P").mkdir(parents=True)
    (tmp_path / "P" / "order.yml").write_text(
        "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
    )
    reference_dir = tmp_path / "P" / "Fam" / "rec1"
    reference_dir.mkdir(parents=True)
    (reference_dir / "proteins.faa").write_text(
        ">rec1|gene_index=0|name=mfa1|role=core_MAT\n" + "M" * 41 + "\n"
        ">rec1|gene_index=1|name=pra1|role=core_MAT\n" + "M" * 300 + "\n"
    )

    def fake_fast_path(proteome_fasta, families, reference_fasta, runner=None):
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    def fake_genomic(*args, **kwargs):
        return []  # mfa1 genuinely not found even after the relaxed second pass

    results = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=fake_genomic,
    )
    assert results[0].genes_missing == []
    assert results[0].genes_not_searchable == ["mfa1"]


def test_short_orf_split_discriminates_three_buckets(tmp_path):
    """Three core genes, three distinct fates: mfa1 is short (41 aa) and
    genuinely never found -> genes_not_searchable. pra2 is normal-length and
    genuinely never found -> genes_missing. pra1 is found -> genes_found.
    A test with only one gene per bucket can't tell "correctly bucketed" from
    "coincidentally correct because nothing else is missing" -- this proves
    the split logic actually discriminates across all three cases at once."""
    three_gene_family = Family(
        FamilyKey("P", "aLocus"), "pattern", None, "^a[0-9]+$",
        [{"name": "mfa1", "role": "core_MAT"},
         {"name": "pra1", "role": "core_MAT"},
         {"name": "pra2", "role": "core_MAT"}],
        [1],
    )
    (tmp_path / "P").mkdir(parents=True)
    (tmp_path / "P" / "order.yml").write_text(
        "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
        "      - {name: pra2, role: core_MAT}\n"
    )
    reference_dir = tmp_path / "P" / "Fam" / "rec1"
    reference_dir.mkdir(parents=True)
    (reference_dir / "proteins.faa").write_text(
        ">rec1|gene_index=0|name=mfa1|role=core_MAT\n" + "M" * 41 + "\n"
        ">rec1|gene_index=1|name=pra1|role=core_MAT\n" + "M" * 300 + "\n"
        ">rec1|gene_index=2|name=pra2|role=core_MAT\n" + "M" * 300 + "\n"
    )

    def fake_fast_path(proteome_fasta, families, reference_fasta, runner=None):
        return [SearchHit(three_gene_family.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1",
                           "diamond_proteome")]

    def fake_genomic(*args, **kwargs):
        return []  # neither mfa1 nor pra2 is found even after the relaxed second pass

    results = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=fake_genomic,
        ambiguity_floor=0.3,  # only 1 of 3 genes is found by design -- lower the floor
        # so this single-family cluster still clears scoring and isn't dropped,
        # without needing a fourth gene just to satisfy an unrelated threshold.
    )
    assert results[0].genes_found == ["pra1"]
    assert results[0].genes_missing == ["pra2"]
    assert results[0].genes_not_searchable == ["mfa1"]
