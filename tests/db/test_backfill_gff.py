from __future__ import annotations

from unittest.mock import MagicMock

import yaml

from MATPredict.db.cli import build_gff_for_record, find_records_missing_proteins_faa


def _write_record(db_root, phylum, order_or_family, record_id, genes=None):
    record_dir = db_root / phylum / order_or_family / record_id
    record_dir.mkdir(parents=True)
    metadata = {
        "record_id": record_id,
        "genes": genes if genes is not None else [
            {"gene_index": 0, "name": "G1", "present": True, "protein_accession": "ncbi_protein:ABC1.1"},
        ],
    }
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(metadata))
    return record_dir


def test_find_records_missing_proteins_faa_skips_records_that_already_have_one(tmp_path):
    db_root = tmp_path / "db"
    _write_record(db_root, "Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")
    complete_dir = _write_record(db_root, "Ascomycota", "Eurotiales", "222_b_MAT_MAT1-2")
    (complete_dir / "proteins.faa").write_text(">already|gene_index=0|name=G1|role=core_MAT\nMSEQ\n")

    missing = find_records_missing_proteins_faa(db_root)

    assert missing == [("Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")]


def test_find_records_missing_proteins_faa_treats_empty_file_as_missing(tmp_path):
    db_root = tmp_path / "db"
    record_dir = _write_record(db_root, "Ascomycota", "Teloschistales", "111_a_MAT_MAT1-1")
    (record_dir / "proteins.faa").write_text("\n")  # empty -- zero FASTA entries

    missing = find_records_missing_proteins_faa(db_root)

    assert missing == [("Ascomycota", "Teloschistales", "111_a_MAT_MAT1-1")]


def test_find_records_missing_proteins_faa_excludes_candidates(tmp_path):
    db_root = tmp_path / "db"
    _write_record(db_root, "candidates", "Ascomycota", "999_c_MAT_MAT1-1")

    assert find_records_missing_proteins_faa(db_root) == []


def test_backfill_missing_proteins_faa_isolates_one_record_failure(tmp_path, monkeypatch):
    db_root = tmp_path / "db"
    _write_record(db_root, "Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")
    _write_record(db_root, "Ascomycota", "Onygenales", "222_b_MAT_MAT1-2")

    calls = []

    def fake_build(db_root, phylum, order_or_family, record_id, ncbi, uniprot):
        calls.append(record_id)
        if record_id == "111_a_MAT_MAT1-1":
            raise ConnectionError("simulated NCBI outage")

    monkeypatch.setattr("MATPredict.db.cli.build_gff_for_record", fake_build)

    from MATPredict.db.cli import backfill_missing_proteins_faa

    succeeded, failed = backfill_missing_proteins_faa(db_root, ncbi=MagicMock(), uniprot=MagicMock())

    assert calls == ["111_a_MAT_MAT1-1", "222_b_MAT_MAT1-2"]  # second record still attempted
    assert succeeded == [("Ascomycota", "Onygenales", "222_b_MAT_MAT1-2")]
    assert len(failed) == 1
    assert failed[0][0] == ("Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")
    assert "simulated NCBI outage" in failed[0][1]


def _write_record_with_coordinates_no_accession(db_root, phylum, order_or_family, record_id):
    """A gene with real segment coordinates but protein_accession: null -- the
    Xanthoria MAG case: a real genomic span exists, no NCBI protein record does."""
    record_dir = db_root / phylum / order_or_family / record_id
    record_dir.mkdir(parents=True)
    metadata = {
        "record_id": record_id,
        "locus": {"core": {"segments": [
            {
                "segment_index": 0,
                "sequence_source": {"type": "insdc_nucleotide", "accession": "ACC1.1", "seq_region": "ACC1.1"},
                "start": 100, "end": 108,
            },
        ]}},
        "genes": [
            {
                "gene_index": 0, "name": "G1", "role": "core_MAT", "present": True,
                "protein_accession": None, "segment_index": 0,
                "start": 100, "end": 108, "strand": "+",
            },
        ],
    }
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(metadata))
    return record_dir


def test_build_gff_for_record_derives_sequence_from_coordinates_when_no_accession(tmp_path, monkeypatch):
    db_root = tmp_path / "db"
    record_dir = _write_record_with_coordinates_no_accession(
        db_root, "Ascomycota", "Teloschistales", "111_a_MAT_combined"
    )

    monkeypatch.setattr(
        "MATPredict.db.cli._independent_translation",
        lambda record, gene, ncbi: "MSEQ",
    )

    build_gff_for_record(db_root, "Ascomycota", "Teloschistales", "111_a_MAT_combined", ncbi=MagicMock(), uniprot=MagicMock())

    faa_text = (record_dir / "proteins.faa").read_text()
    assert ">111_a_MAT_combined|gene_index=0|name=G1|role=core_MAT" in faa_text
    assert "MSEQ" in faa_text


def test_build_gff_for_record_skips_gene_when_translation_returns_none(tmp_path, monkeypatch):
    db_root = tmp_path / "db"
    _write_record_with_coordinates_no_accession(db_root, "Ascomycota", "Teloschistales", "222_b_MAT_combined")

    monkeypatch.setattr("MATPredict.db.cli._independent_translation", lambda record, gene, ncbi: None)

    build_gff_for_record(db_root, "Ascomycota", "Teloschistales", "222_b_MAT_combined", ncbi=MagicMock(), uniprot=MagicMock())

    faa_text = (db_root / "Ascomycota" / "Teloschistales" / "222_b_MAT_combined" / "proteins.faa").read_text()
    assert faa_text.strip() == ""


# --- staleness selector: find records whose locus.gbk predates the real-CDS generator ---

def _write_locus_gbk(record_dir, feature_types):
    """Write a minimal but genuinely parseable locus.gbk carrying one feature per
    entry in `feature_types` ("gene", "CDS", ...). Written with Bio.SeqIO so the
    file's column layout is whatever Biopython really produces, rather than a
    hand-built string the parser might disagree with."""
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqIO import write as seqio_write
    from Bio.SeqRecord import SeqRecord

    seq_record = SeqRecord(Seq("ATG" * 10), id="seg0", name="seg0", description="test segment")
    seq_record.annotations["molecule_type"] = "DNA"
    for feature_type in feature_types:
        qualifiers = {"gene": ["G1"]}
        if feature_type == "CDS":
            qualifiers["translation"] = ["MSEQ"]
        seq_record.features.append(
            SeqFeature(FeatureLocation(0, 30, strand=1), type=feature_type, qualifiers=qualifiers)
        )
    seqio_write([seq_record], record_dir / "locus.gbk", "genbank")


def test_find_records_with_stale_locus_gbk_skips_record_with_real_cds_features(tmp_path):
    from MATPredict.db.cli import find_records_with_stale_locus_gbk

    db_root = tmp_path / "db"
    record_dir = _write_record(db_root, "Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")
    _write_locus_gbk(record_dir, ["gene", "CDS"])

    assert find_records_with_stale_locus_gbk(db_root) == []


def test_find_records_with_stale_locus_gbk_returns_record_with_only_gene_features(tmp_path):
    """The pre-upgrade generator wrote `gene` features and an all-N ORIGIN block but
    no `CDS`/`translation` -- that is exactly the staleness this selector must catch."""
    from MATPredict.db.cli import find_records_with_stale_locus_gbk

    db_root = tmp_path / "db"
    record_dir = _write_record(db_root, "Ascomycota", "Eurotiales", "222_b_MAT_MAT1-2")
    _write_locus_gbk(record_dir, ["gene"])

    assert find_records_with_stale_locus_gbk(db_root) == [
        ("Ascomycota", "Eurotiales", "222_b_MAT_MAT1-2")
    ]


def test_find_records_with_stale_locus_gbk_returns_record_with_no_locus_gbk(tmp_path):
    from MATPredict.db.cli import find_records_with_stale_locus_gbk

    db_root = tmp_path / "db"
    _write_record(db_root, "Ascomycota", "Onygenales", "333_c_MAT_MAT1-1")

    assert find_records_with_stale_locus_gbk(db_root) == [
        ("Ascomycota", "Onygenales", "333_c_MAT_MAT1-1")
    ]


def test_find_records_with_stale_locus_gbk_excludes_candidates(tmp_path):
    from MATPredict.db.cli import find_records_with_stale_locus_gbk

    db_root = tmp_path / "db"
    _write_record(db_root, "candidates", "Ascomycota", "999_c_MAT_MAT1-1")

    assert find_records_with_stale_locus_gbk(db_root) == []


def test_find_records_with_stale_locus_gbk_skips_record_whose_genes_are_all_absent(tmp_path):
    """A record whose every gene is `present: false` has nothing to regenerate: a
    freshly built locus.gbk would still carry zero CDS features, so returning it
    would make the sweep loop on it forever."""
    from MATPredict.db.cli import find_records_with_stale_locus_gbk

    db_root = tmp_path / "db"
    _write_record(
        db_root, "Ascomycota", "Eurotiales", "444_d_MAT_MAT1-2",
        genes=[{"gene_index": 0, "name": "G1", "present": False, "protein_accession": None}],
    )

    assert find_records_with_stale_locus_gbk(db_root) == []


def test_find_records_with_stale_locus_gbk_treats_unparseable_file_as_stale(tmp_path):
    """A locus.gbk Bio.SeqIO cannot parse is not evidence of freshness; regenerating
    it is the fix, so it is selected rather than crashing the sweep."""
    from MATPredict.db.cli import find_records_with_stale_locus_gbk

    db_root = tmp_path / "db"
    record_dir = _write_record(db_root, "Ascomycota", "Eurotiales", "555_e_MAT_MAT1-1")
    (record_dir / "locus.gbk").write_text("this is not a GenBank file at all\n")

    assert find_records_with_stale_locus_gbk(db_root) == [
        ("Ascomycota", "Eurotiales", "555_e_MAT_MAT1-1")
    ]


def test_backfill_stale_locus_gbk_isolates_one_record_failure(tmp_path, monkeypatch):
    """The stale-gbk selector drives the same failure-isolating loop the
    proteins.faa selector already does -- one record's live-fetch error never
    aborts the rest of the sweep."""
    from MATPredict.db.cli import _backfill_records

    db_root = tmp_path / "db"
    _write_record(db_root, "Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")
    _write_record(db_root, "Ascomycota", "Onygenales", "222_b_MAT_MAT1-2")

    def fake_build(db_root, phylum, order_or_family, record_id, ncbi, uniprot):
        if record_id == "111_a_MAT_MAT1-1":
            raise ConnectionError("simulated NCBI outage")

    monkeypatch.setattr("MATPredict.db.cli.build_gff_for_record", fake_build)

    from MATPredict.db.cli import find_records_with_stale_locus_gbk

    succeeded, failed = _backfill_records(
        db_root, find_records_with_stale_locus_gbk(db_root), ncbi=MagicMock(), uniprot=MagicMock()
    )

    assert succeeded == [("Ascomycota", "Onygenales", "222_b_MAT_MAT1-2")]
    assert len(failed) == 1
    assert failed[0][0] == ("Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")


# --- placeholder reporting: name the segments that still lack real nucleotide sequence ---

def _write_locus_gbk_via_write_genbank(record_dir, ncbi=None):
    """Write a locus.gbk through the real gff_export.write_genbank, so the file under
    test is exactly what the backfill itself produces -- including how it names each
    segment SeqRecord -- rather than a hand-built stand-in."""
    from MATPredict.db.gff_export import write_genbank

    record = {
        "record_id": record_dir.name,
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 1, "end": 30,
             "sequence_source": {"type": "insdc_nucleotide", "accession": "ACC1.1",
                                 "seq_region": "ACC1.1"}},
        ]}},
        "genes": [],
    }
    write_genbank(record, sequences={}, out_path=record_dir / "locus.gbk", ncbi=ncbi)
    return record_dir / "locus.gbk"


def test_placeholder_segments_names_an_all_n_segment(tmp_path):
    """A segment whose nucleotide fetch could not happen is written as all-N. The
    sweep must name it, so the operator learns from the run's own output which
    records still lack real sequence."""
    from MATPredict.db.cli import placeholder_segments

    db_root = tmp_path / "db"
    record_dir = _write_record(db_root, "Ascomycota", "Eurotiales", "666_f_MAT_MAT1-1")
    gbk_path = _write_locus_gbk_via_write_genbank(record_dir)  # no ncbi -- all-N placeholder

    assert placeholder_segments(gbk_path) == ["666_f_MAT_MAT1-1.segment0"]


def test_placeholder_segments_returns_empty_for_real_sequence(tmp_path):
    from MATPredict.db.cli import placeholder_segments

    db_root = tmp_path / "db"
    record_dir = _write_record(db_root, "Ascomycota", "Eurotiales", "777_g_MAT_MAT1-2")
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.return_value = "ACGT" * 7 + "TA"
    gbk_path = _write_locus_gbk_via_write_genbank(record_dir, ncbi=fake_ncbi)

    assert placeholder_segments(gbk_path) == []


def test_placeholder_segments_returns_empty_for_an_unreadable_file(tmp_path):
    """Reporting must never turn a succeeded record into a failure."""
    from MATPredict.db.cli import placeholder_segments

    gbk_path = tmp_path / "locus.gbk"
    gbk_path.write_text("this is not a GenBank file at all\n")

    assert placeholder_segments(gbk_path) == []
