from __future__ import annotations

import yaml

from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.detect.benchmark import (
    match_ground_truth,
    run_benchmark,
    score_self_consistency,
)


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


# --- Task 5: self-consistency ground-truth scoring -------------------------

# A curated record with one present gene ("g1"), sourced from accession
# TESTACC.1 -- deliberately independent of the pilot rollout's real
# Coccidioides/Aspergillus taxids/accessions so this fixture never
# accidentally collides with real db/ content.
_CURATED_TAXID = 5501
_CURATED_RECORD_ID = "5501_teststrain_L_a1"


def _write_curated_db(tmp_path, source_accession="TESTACC.1", strain="teststrain"):
    (tmp_path / "P").mkdir(exist_ok=True)
    order_yml = tmp_path / "P" / "order.yml"
    if not order_yml.exists():
        order_yml.write_text(
            "phylum: P\nloci:\n  - locus_name: L\n    vocabulary_type: pattern\n"
            "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
            "    genes:\n      - {name: g1, role: core_MAT}\n"
        )
    record_dir = tmp_path / "P" / "Fam" / _CURATED_RECORD_ID
    record_dir.mkdir(parents=True, exist_ok=True)
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump({
        "record_id": _CURATED_RECORD_ID,
        "taxonomy": {"taxid": _CURATED_TAXID, "lineage": "k__Fungi;g__X;s__X_sp"},
        "organism": {"species": "X sp", "strain": {"name": strain}},
        "mating_type": {"locus_name": "L", "idiomorphs": ["a1"]},
        "locus": {
            "coordinate_provenance": "published_explicit",
            "excluded_from_coordinate_benchmark": False,
            "core": {"segments": [{
                "sequence_source": {
                    "type": "insdc_nucleotide", "accession": source_accession,
                    "seq_region": source_accession,
                },
                "start": 1, "end": 100,
            }]},
        },
        "genes": [{
            "gene_index": 0, "name": "g1", "role": "core_MAT", "present": True,
            "segment_index": 0, "start": 1, "end": 50, "strand": "+",
            "protein_accession": "ncbi_protein:CURATED1.1",
        }],
    }, sort_keys=False))
    return record_dir


def _write_rollout_report(tmp_path, genome_id, contig, start, end, strand):
    genome_dir = tmp_path / "rollout" / genome_id
    genome_dir.mkdir(parents=True, exist_ok=True)
    report_path = genome_dir / "detection_report.yaml"
    report_path.write_text(yaml.safe_dump({
        "families_attempted": ["P:L"],
        "detected": [{
            "family": "P:L",
            "gene_evidence": [{
                "gene": "g1", "role": "core_MAT", "contig": contig,
                "start": start, "end": end, "strand": strand,
                "identity": 100.0, "coverage": 100.0,
                "reference_record": None, "method": "exonerate", "status": "pass",
                "alternate_model": None,
            }],
        }],
        "not_detected": [],
    }, sort_keys=False))
    return report_path


def _write_genome_fasta(tmp_path, name, contig, sequence):
    fasta_path = tmp_path / name
    fasta_path.write_text(f">{contig}\n{sequence}\n")
    return fasta_path


def _fake_ncbi(protein_fasta_by_accession, tmp_path):
    def transport(url: str) -> str:
        for accession, body in protein_fasta_by_accession.items():
            if f"id={accession}" in url:
                return body
        raise AssertionError(f"unexpected url: {url}")
    fetcher = CachedFetcher(cache_dir=tmp_path / "ncbi_cache", transport=transport)
    return NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)


# ATGAAATAA -> translate(to_stop=True) -> "MK" (stop codon TAA dropped).
_ROLLOUT_NUC = "ATGAAATAA"
_ROLLOUT_PROTEIN_FASTA = ">CURATED1.1\nMK\n"
_DIFFERENT_PROTEIN_FASTA = ">CURATED1.1\nQWERTYASDFGH\n"


def test_score_self_consistency_exact_protein_match_scores_as_found(tmp_path):
    _write_curated_db(tmp_path)
    genome_id = f"{_CURATED_TAXID}_TESTACC.1"
    report_path = _write_rollout_report(tmp_path, genome_id, "ctg1", 1, 9, "+")
    fasta_path = _write_genome_fasta(tmp_path, "genome.fa", "ctg1", _ROLLOUT_NUC)
    ncbi = _fake_ncbi({"CURATED1.1": _ROLLOUT_PROTEIN_FASTA}, tmp_path)

    scored, matches = score_self_consistency(
        report_paths=[report_path], db_root=tmp_path,
        genome_fasta_paths={genome_id: fasta_path}, ncbi=ncbi,
    )

    assert len(scored) == 1
    assert scored[0].sensitivity == 1.0
    assert len(matches) == 1
    assert matches[0].status == "exact"


def test_score_self_consistency_wrong_protein_scores_as_not_found(tmp_path):
    _write_curated_db(tmp_path)
    genome_id = f"{_CURATED_TAXID}_TESTACC.1"
    report_path = _write_rollout_report(tmp_path, genome_id, "ctg1", 1, 9, "+")
    fasta_path = _write_genome_fasta(tmp_path, "genome.fa", "ctg1", _ROLLOUT_NUC)
    ncbi = _fake_ncbi({"CURATED1.1": _DIFFERENT_PROTEIN_FASTA}, tmp_path)

    scored, matches = score_self_consistency(
        report_paths=[report_path], db_root=tmp_path,
        genome_fasta_paths={genome_id: fasta_path}, ncbi=ncbi,
    )

    assert len(scored) == 1
    assert scored[0].sensitivity == 0.0
    assert matches[0].status == "exact"


def test_score_self_consistency_ambiguous_strain_excluded_and_flagged(tmp_path):
    _write_curated_db(tmp_path, source_accession="TESTACC.1", strain="teststrain")
    # A rollout genome with the SAME taxid but a DIFFERENT accession than the
    # curated record's own source -- the real-world case this project hit for
    # every one of its pilot Coccidioides/Aspergillus genomes.
    genome_id = f"{_CURATED_TAXID}_OTHERACC.1"
    report_path = _write_rollout_report(tmp_path, genome_id, "ctg1", 1, 9, "+")
    fasta_path = _write_genome_fasta(tmp_path, "genome.fa", "ctg1", _ROLLOUT_NUC)
    ncbi = _fake_ncbi({"CURATED1.1": _ROLLOUT_PROTEIN_FASTA}, tmp_path)

    scored, matches = score_self_consistency(
        report_paths=[report_path], db_root=tmp_path,
        genome_fasta_paths={genome_id: fasta_path}, ncbi=ncbi,
    )

    # Not scored numerically...
    assert scored == []
    # ...but flagged, not silently dropped.
    assert len(matches) == 1
    assert matches[0].status == "ambiguous"
    assert "same taxid" in matches[0].reason


def test_match_ground_truth_no_curated_record_for_taxid_returns_empty(tmp_path):
    _write_curated_db(tmp_path)
    assert match_ground_truth("999999_SOMEACC.1", tmp_path) == []
