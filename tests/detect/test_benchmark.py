from __future__ import annotations

import yaml

from MATPredict.db import validate as db_validate
from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.detect.benchmark import (
    _extract_translated_gene,
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
_CURATED_TAXID_B = 199306
_CURATED_RECORD_ID_B = "199306_teststrain_L_a1"

_DEFAULT_GENES = [{
    "gene_index": 0, "name": "g1", "role": "core_MAT", "present": True,
    "segment_index": 0, "start": 1, "end": 50, "strand": "+",
    "protein_accession": "ncbi_protein:CURATED1.1",
}]


def _missing_manifest(tmp_path):
    """A manifest path that never exists on disk, for full isolation from
    the real BFD acquisition manifest mounted on UCR HPCC
    (`genome_acquisition.DEFAULT_LOCAL_MANIFEST_PATH`) -- deterministic
    regardless of what environment the suite runs in."""
    return tmp_path / "no_such_manifest.csv"


def _write_manifest(tmp_path, rows):
    """A minimal BFD-manifest-shaped CSV (same columns `_manifest_strain`
    reads: ASMID, STRAIN, NCBI_TAXONID) for strain-match tests."""
    manifest_path = tmp_path / "manifest.csv"
    header = "ASMID,SPECIES_IN,STRAIN,NCBI_TAXONID\n"
    body = "".join(f"{r['ASMID']},{r.get('SPECIES_IN','')},{r['STRAIN']},{r['NCBI_TAXONID']}\n" for r in rows)
    manifest_path.write_text(header + body)
    return manifest_path


def _write_curated_db(
    tmp_path, taxid=_CURATED_TAXID, record_id=_CURATED_RECORD_ID,
    source_accession="TESTACC.1", strain="teststrain", genes=None,
):
    (tmp_path / "P").mkdir(exist_ok=True)
    order_yml = tmp_path / "P" / "order.yml"
    if not order_yml.exists():
        order_yml.write_text(
            "phylum: P\nloci:\n  - locus_name: L\n    vocabulary_type: pattern\n"
            "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
            "    genes:\n      - {name: g1, role: core_MAT}\n"
        )
    record_dir = tmp_path / "P" / "Fam" / record_id
    record_dir.mkdir(parents=True, exist_ok=True)
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump({
        "record_id": record_id,
        "taxonomy": {"taxid": taxid, "lineage": "k__Fungi;g__X;s__X_sp"},
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
        "genes": genes if genes is not None else _DEFAULT_GENES,
    }, sort_keys=False))
    return record_dir


def _write_rollout_report(tmp_path, genome_id, gene_evidence):
    """`gene_evidence` is a list of dicts, each at least
    `{gene, contig, start, end, strand}`; `exons` (a list of
    `{start, end}` dicts) and `status`/`method`/`identity`/`coverage` are
    optional, matching `report.py`'s real `_result_doc` shape."""
    genome_dir = tmp_path / "rollout" / genome_id
    genome_dir.mkdir(parents=True, exist_ok=True)
    report_path = genome_dir / "detection_report.yaml"
    entries = [{
        "gene": ev["gene"], "role": ev.get("role", "core_MAT"), "contig": ev["contig"],
        "start": ev["start"], "end": ev["end"], "strand": ev.get("strand"),
        "identity": ev.get("identity", 100.0), "coverage": ev.get("coverage", 100.0),
        "reference_record": None, "method": ev.get("method", "exonerate"),
        "status": ev.get("status", "pass"), "alternate_model": None,
        "exons": ev.get("exons"),
    } for ev in gene_evidence]
    report_path.write_text(yaml.safe_dump({
        "families_attempted": ["P:L"],
        "detected": [{"family": "P:L", "gene_evidence": entries}],
        "not_detected": [],
    }, sort_keys=False))
    return report_path


def _write_genome_fasta(tmp_path, name, contig, sequence):
    fasta_path = tmp_path / name
    fasta_path.write_text(f">{contig}\n{sequence}\n")
    return fasta_path


def _fake_ncbi(protein_fasta_by_accession, tmp_path, cache_name="ncbi_cache"):
    def transport(url: str) -> str:
        for accession, body in protein_fasta_by_accession.items():
            if f"id={accession}" in url:
                return body
        raise AssertionError(f"unexpected url: {url}")
    fetcher = CachedFetcher(cache_dir=tmp_path / cache_name, transport=transport)
    return NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)


# ATGAAATAA -> translate(to_stop=True) -> "MK" (stop codon TAA dropped).
_ROLLOUT_NUC = "ATGAAATAA"
_ROLLOUT_PROTEIN_FASTA = ">CURATED1.1\nMK\n"
_DIFFERENT_PROTEIN_FASTA = ">CURATED1.1\nQWERTYASDFGH\n"


def test_score_self_consistency_exact_protein_match_scores_as_found(tmp_path):
    _write_curated_db(tmp_path)
    genome_id = f"{_CURATED_TAXID}_TESTACC.1"
    report_path = _write_rollout_report(
        tmp_path, genome_id, [{"gene": "g1", "contig": "ctg1", "start": 1, "end": 9, "strand": "+"}],
    )
    fasta_path = _write_genome_fasta(tmp_path, "genome.fa", "ctg1", _ROLLOUT_NUC)
    ncbi = _fake_ncbi({"CURATED1.1": _ROLLOUT_PROTEIN_FASTA}, tmp_path)

    scored, matches = score_self_consistency(
        report_paths=[report_path], db_root=tmp_path,
        genome_fasta_paths={genome_id: fasta_path}, ncbi=ncbi,
        manifest_path=_missing_manifest(tmp_path),
    )

    assert len(scored) == 1
    assert scored[0].sensitivity == 1.0
    assert len(matches) == 1
    assert matches[0].status == "exact"


def test_score_self_consistency_wrong_protein_scores_as_not_found(tmp_path):
    _write_curated_db(tmp_path)
    genome_id = f"{_CURATED_TAXID}_TESTACC.1"
    report_path = _write_rollout_report(
        tmp_path, genome_id, [{"gene": "g1", "contig": "ctg1", "start": 1, "end": 9, "strand": "+"}],
    )
    fasta_path = _write_genome_fasta(tmp_path, "genome.fa", "ctg1", _ROLLOUT_NUC)
    ncbi = _fake_ncbi({"CURATED1.1": _DIFFERENT_PROTEIN_FASTA}, tmp_path)

    scored, matches = score_self_consistency(
        report_paths=[report_path], db_root=tmp_path,
        genome_fasta_paths={genome_id: fasta_path}, ncbi=ncbi,
        manifest_path=_missing_manifest(tmp_path),
    )

    assert len(scored) == 1
    assert scored[0].sensitivity == 0.0
    assert matches[0].status == "exact"


def test_score_self_consistency_ambiguous_strain_excluded_and_flagged(tmp_path):
    _write_curated_db(tmp_path, source_accession="TESTACC.1", strain="teststrain")
    # A rollout genome with the SAME taxid but a DIFFERENT accession than the
    # curated record's own source, and no manifest strain info available --
    # the real-world case this project hit for every one of its pilot
    # Coccidioides/Aspergillus genomes.
    genome_id = f"{_CURATED_TAXID}_OTHERACC.1"
    report_path = _write_rollout_report(
        tmp_path, genome_id, [{"gene": "g1", "contig": "ctg1", "start": 1, "end": 9, "strand": "+"}],
    )
    fasta_path = _write_genome_fasta(tmp_path, "genome.fa", "ctg1", _ROLLOUT_NUC)
    ncbi = _fake_ncbi({"CURATED1.1": _ROLLOUT_PROTEIN_FASTA}, tmp_path)

    scored, matches = score_self_consistency(
        report_paths=[report_path], db_root=tmp_path,
        genome_fasta_paths={genome_id: fasta_path}, ncbi=ncbi,
        manifest_path=_missing_manifest(tmp_path),
    )

    # Not scored numerically...
    assert scored == []
    # ...but flagged, not silently dropped.
    assert len(matches) == 1
    assert matches[0].status == "ambiguous"
    assert "same taxid" in matches[0].reason


def test_match_ground_truth_no_curated_record_for_taxid_returns_empty(tmp_path):
    _write_curated_db(tmp_path)
    assert match_ground_truth(
        "999999_SOMEACC.1", tmp_path, manifest_path=_missing_manifest(tmp_path),
    ) == []


# --- Fix-round: strain-aware matching ---------------------------------------

def test_match_ground_truth_strain_match_scores_as_exact(tmp_path):
    """Same taxid, DIFFERENT accession, but the BFD manifest's own STRAIN
    column for this rollout genome equals the curated record's strain --
    the real 'same strain, different assembly' case the brief's Step 2
    asked for (e.g. a genuine FGSC A4 rollout genome vs. the curated FGSC A4
    A. nidulans record) that pure accession matching alone would wrongly
    mark ambiguous."""
    _write_curated_db(tmp_path, source_accession="TESTACC.1", strain="H538.4")
    manifest_path = _write_manifest(tmp_path, [{
        "ASMID": "OTHERACC.1_some_assembly_name", "STRAIN": "H538.4",
        "NCBI_TAXONID": str(_CURATED_TAXID),
    }])

    matches = match_ground_truth(
        f"{_CURATED_TAXID}_OTHERACC.1", tmp_path, manifest_path=manifest_path,
    )

    assert len(matches) == 1
    assert matches[0].status == "exact"
    assert "strain" in matches[0].reason.lower()


def test_match_ground_truth_strain_mismatch_stays_ambiguous(tmp_path):
    """Manifest strain IS available but genuinely differs from the curated
    record's strain -- must stay ambiguous, not be waved through just
    because a strain lookup succeeded at all."""
    _write_curated_db(tmp_path, source_accession="TESTACC.1", strain="H538.4")
    manifest_path = _write_manifest(tmp_path, [{
        "ASMID": "OTHERACC.1_some_assembly_name", "STRAIN": "WA_211",
        "NCBI_TAXONID": str(_CURATED_TAXID),
    }])

    matches = match_ground_truth(
        f"{_CURATED_TAXID}_OTHERACC.1", tmp_path, manifest_path=manifest_path,
    )

    assert len(matches) == 1
    assert matches[0].status == "ambiguous"


def test_match_ground_truth_accession_match_scores_exact_despite_taxid_mismatch(tmp_path):
    """BFD carries strain-level taxids for some rollout genomes (e.g.
    Agaricus bisporus var. bisporus = 192523) that differ from a curated
    record's species-level taxid (e.g. A. bisporus = 5346), even though the
    rollout genome's own accession is literally the curated record's own
    source accession. The taxid gate used to run BEFORE the accession check,
    so this genome -- unambiguous ground truth by accession -- was silently
    skipped entirely (matched against zero records), the real defect found
    scoring the Basidiomycota order-testing rollout (3/10 curated records
    reachable; 2 recovered once this is fixed). Accession identity is
    unconditional evidence regardless of which taxid either side declares."""
    _write_curated_db(tmp_path, taxid=5346, source_accession="GCF_000143185.2")

    matches = match_ground_truth(
        "192523_GCF_000143185.2", tmp_path, manifest_path=_missing_manifest(tmp_path),
    )

    assert len(matches) == 1
    assert matches[0].status == "exact"
    assert matches[0].record_id == _CURATED_RECORD_ID
    assert "accession" in matches[0].reason.lower()


def test_match_ground_truth_taxid_mismatch_with_no_accession_match_returns_nothing(tmp_path):
    """The flip side: a taxid mismatch with no accession evidence either must
    still return no match -- the fix must not turn into 'ignore taxid
    entirely' and start reporting unrelated species as ground truth."""
    _write_curated_db(tmp_path, taxid=5346, source_accession="GCF_000143185.2")

    matches = match_ground_truth(
        "999999_UNRELATEDACC.1", tmp_path, manifest_path=_missing_manifest(tmp_path),
    )

    assert matches == []


# --- Fix-round: not-evaluable genes excluded, never fabricated as misses ---

_TWO_GENES = [
    {
        "gene_index": 0, "name": "g1", "role": "core_MAT", "present": True,
        "segment_index": 0, "start": 1, "end": 50, "strand": "+",
        "protein_accession": "ncbi_protein:CURATED1.1",
    },
    {
        "gene_index": 1, "name": "g2", "role": "core_MAT", "present": True,
        "segment_index": 0, "start": 60, "end": 100, "strand": "+",
        "protein_accession": "ncbi_protein:CURATED2.1",
    },
]


def test_score_self_consistency_missing_contig_excluded_not_fabricated(tmp_path):
    """g1 has real gene_evidence and a matching genome FASTA contig (found).
    g2 has gene_evidence too, but its reported contig doesn't exist in the
    supplied genome FASTA (a real-world case: contig-naming mismatch)
    -- that must be excluded from the denominator, not scored as a miss
    that would drag sensitivity down to 0.5."""
    _write_curated_db(tmp_path, genes=_TWO_GENES)
    genome_id = f"{_CURATED_TAXID}_TESTACC.1"
    report_path = _write_rollout_report(tmp_path, genome_id, [
        {"gene": "g1", "contig": "ctg1", "start": 1, "end": 9, "strand": "+"},
        {"gene": "g2", "contig": "NOT_A_REAL_CONTIG", "start": 1, "end": 9, "strand": "+"},
    ])
    fasta_path = _write_genome_fasta(tmp_path, "genome.fa", "ctg1", _ROLLOUT_NUC)
    ncbi = _fake_ncbi({"CURATED1.1": _ROLLOUT_PROTEIN_FASTA}, tmp_path)

    scored, _matches = score_self_consistency(
        report_paths=[report_path], db_root=tmp_path,
        genome_fasta_paths={genome_id: fasta_path}, ncbi=ncbi,
        manifest_path=_missing_manifest(tmp_path),
    )

    assert len(scored) == 1
    assert scored[0].sensitivity == 1.0  # g2 excluded, not counted as a miss
    assert "not evaluable" in scored[0].note
    assert "g2" in scored[0].note


def test_score_self_consistency_missing_genome_fasta_excludes_without_fabricating(tmp_path):
    """No genome FASTA at all is supplied for this genome -- must not
    fabricate a sensitivity=0.0 attributed to the pipeline; the pairing is
    still an "exact" ground-truth match, just not evaluable, so nothing is
    scored for it at all."""
    _write_curated_db(tmp_path)
    genome_id = f"{_CURATED_TAXID}_TESTACC.1"
    report_path = _write_rollout_report(
        tmp_path, genome_id, [{"gene": "g1", "contig": "ctg1", "start": 1, "end": 9, "strand": "+"}],
    )
    ncbi = _fake_ncbi({"CURATED1.1": _ROLLOUT_PROTEIN_FASTA}, tmp_path)

    scored, matches = score_self_consistency(
        report_paths=[report_path], db_root=tmp_path,
        genome_fasta_paths={},  # no FASTA supplied for this genome
        ncbi=ncbi, manifest_path=_missing_manifest(tmp_path),
    )

    assert scored == []
    assert matches[0].status == "exact"  # still a real match, just unscorable


# --- Fix-round: a live NCBI fetch failure is excluded, not fabricated, and
# never aborts the rest of the batch --------------------------------------

def test_score_self_consistency_ncbi_fetch_failure_excluded_and_does_not_abort_batch(tmp_path):
    _write_curated_db(tmp_path, taxid=_CURATED_TAXID, record_id=_CURATED_RECORD_ID)
    _write_curated_db(
        tmp_path, taxid=_CURATED_TAXID_B, record_id=_CURATED_RECORD_ID_B,
        source_accession="OKACC.1",
        genes=[{
            "gene_index": 0, "name": "g1", "role": "core_MAT", "present": True,
            "segment_index": 0, "start": 1, "end": 50, "strand": "+",
            "protein_accession": "ncbi_protein:CURATED2.1",
        }],
    )

    # Genome A's protein accession (CURATED1.1) is deliberately NOT in the
    # fake NCBI transport's map, so fetch_protein_sequence raises for it.
    genome_id_a = f"{_CURATED_TAXID}_TESTACC.1"
    report_path_a = _write_rollout_report(
        tmp_path, genome_id_a, [{"gene": "g1", "contig": "ctg1", "start": 1, "end": 9, "strand": "+"}],
    )
    fasta_path_a = _write_genome_fasta(tmp_path, "genome_a.fa", "ctg1", _ROLLOUT_NUC)

    # Genome B's own protein accession (CURATED2.1) IS in the map and
    # matches -- this genome must still be scored even though genome A's
    # fetch failed.
    genome_id_b = f"{_CURATED_TAXID_B}_OKACC.1"
    report_path_b = _write_rollout_report(
        tmp_path, genome_id_b, [{"gene": "g1", "contig": "ctg2", "start": 1, "end": 9, "strand": "+"}],
    )
    fasta_path_b = _write_genome_fasta(tmp_path, "genome_b.fa", "ctg2", _ROLLOUT_NUC)

    ncbi = _fake_ncbi({"CURATED2.1": _ROLLOUT_PROTEIN_FASTA}, tmp_path)

    scored, matches = score_self_consistency(
        report_paths=[report_path_a, report_path_b], db_root=tmp_path,
        genome_fasta_paths={genome_id_a: fasta_path_a, genome_id_b: fasta_path_b},
        ncbi=ncbi, manifest_path=_missing_manifest(tmp_path),
    )

    # Genome A: NCBI fetch failed -> not evaluable -> no fabricated score.
    assert not any(genome_id_a in fb.note for fb in scored)
    genome_a_matches = [m for m in matches if m.genome_id == genome_id_a]
    assert genome_a_matches and genome_a_matches[0].status == "exact"

    # Genome B: unaffected by genome A's failure, scored normally.
    genome_b_scores = [fb for fb in scored if genome_id_b in fb.note]
    assert len(genome_b_scores) == 1
    assert genome_b_scores[0].sensitivity == 1.0

    # Exactly one genome ended up scored (A did not), proving the failure
    # was excluded rather than fabricated as a 0.0 AND that it did not
    # abort processing of genome B.
    assert len(scored) == 1


# --- Fix-round (Critical): multi-exon splicing, cross-checked against
# db/validate.py's own independent transcript-assembly logic ---------------

def test_extract_translated_gene_multi_exon_matches_validate_transcript_assembly(tmp_path):
    """Real 2-exon MAT1-1-1 coordinates (plus strand) from
    db/Ascomycota/Onygenales/199306_rmscc1040_MAT_MAT1-1/metadata.yaml:
    exons [(6905, 7196), (7250, 8127)]. Translating the RAW 6905-8127 span
    (1223 nt, including the 53nt intron) would read through the intron and
    corrupt the protein for a perfectly correct detection -- exactly the
    Critical bug this fix addresses. This test builds a synthetic genome
    whose only real content is an in-frame ORF spliced across exactly these
    two exon spans, and cross-checks `_extract_translated_gene`'s output
    against `db/validate.py`'s own, independently-written transcript
    assembly + translation logic (`_independent_translation`) fed the exact
    same underlying sequence via a mocked NcbiClient -- not a hand-computed
    expected string.
    """
    exons = [(6905, 7196), (7250, 8127)]
    total_len = sum(end - start + 1 for start, end in exons)
    assert total_len % 3 == 0
    # A deterministic in-frame ORF (no stop codons) of exactly the spliced length.
    cds = "ATG" + "AAA" * ((total_len - 3) // 3)
    assert len(cds) == total_len

    contig_len = exons[-1][1]
    contig_chars = ["A"] * contig_len
    cursor = 0
    for start, end in exons:
        length = end - start + 1
        contig_chars[start - 1:end] = list(cds[cursor:cursor + length])
        cursor += length
    contig_seq = "".join(contig_chars)
    fasta_path = _write_genome_fasta(tmp_path, "multi_exon_genome.fa", "ctg1", contig_seq)

    actual = _extract_translated_gene(
        fasta_path, "ctg1", exons[0][0], exons[-1][1], "+", exons,
    )

    # Cross-check: db/validate.py's own transcript-assembly + translation,
    # fed the identical sequence via a mocked NcbiClient (same fake-transport
    # style already used by tests/db/test_validate.py and test_ncbi_client.py).
    def transport(url: str) -> str:
        query = dict(p.split("=", 1) for p in url.split("?", 1)[1].split("&") if "=" in p)
        start, stop = int(query["seq_start"]), int(query["seq_stop"])
        return f">sub\n{contig_seq[start - 1:stop]}\n"

    fetcher = CachedFetcher(cache_dir=tmp_path / "ncbi_cache_validate", transport=transport)
    ncbi = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    record = {"locus": {"core": {"segments": [{
        "sequence_source": {"type": "insdc_nucleotide", "accession": "FAKEACC.1"},
    }]}}}
    gene = {
        "segment_index": 0, "strand": "+", "codon_start": 1, "transl_table": 1,
        "exons": [{"start": s, "end": e} for s, e in exons],
    }
    expected = db_validate._independent_translation(record, gene, ncbi)

    assert expected is not None
    assert actual == expected
    # And sanity: it must NOT equal what naive whole-span (intron-included)
    # translation would produce -- proving the fix actually changed behavior.
    naive = _extract_translated_gene(fasta_path, "ctg1", exons[0][0], exons[-1][1], "+", None)
    assert naive != actual
