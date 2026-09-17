from __future__ import annotations

from MATPredict.detect.benchmark import run_benchmark


def test_run_benchmark_reports_na_for_thin_families(tmp_path):
    (tmp_path / "P").mkdir()
    (tmp_path / "P" / "order.yml").write_text(
        "phylum: P\nloci:\n  - locus_name: L\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: g1, role: core_MAT}\n"
    )
    record_dir = tmp_path / "P" / "Fam" / "1_strain_L_a1"
    record_dir.mkdir(parents=True)
    (record_dir / "metadata.yaml").write_text(
        "record_id: 1_strain_L_a1\n"
        "taxonomy: {taxid: 1, lineage: 'k__Fungi;g__X;s__X_sp'}\n"
        "organism: {species: 'X sp'}\n"
        "mating_type: {locus_name: L, idiomorphs: [a1]}\n"
        "locus: {coordinate_provenance: published_explicit, excluded_from_coordinate_benchmark: false, "
        "core: {segments: [{sequence_source: {type: insdc_nucleotide, accession: 'X.1', seq_region: 'X.1'}, "
        "start: 1, end: 100}]}}\n"
        "genes: [{gene_index: 0, name: g1, role: core_MAT, present: true, segment_index: 0, start: 1, end: 50, strand: '+'}]\n"
    )
    results = run_benchmark(tmp_path)
    assert len(results) == 1
    assert results[0].sensitivity is None
    assert "insufficient data" in results[0].note
    assert results[0].n_reference_after_holdout == 0
