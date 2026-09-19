from __future__ import annotations

from Bio import SeqIO

from MATPredict.db.gff_export import write_genbank, write_gff3, write_proteins_fasta

RECORD = {
    "record_id": "4837_nrrl-1555_MAT_Plus",
    "locus": {
        "core": {
            "segments": [{"segment_index": 0, "sequence_source": {"seq_region": "scaffold_3"}, "start": 120345, "end": 128900}],
        },
    },
    "genes": [
        {"gene_index": 0, "name": "sexP", "role": "core_MAT", "present": True, "segment_index": 0,
         "start": 121002, "end": 122400, "strand": "+"},
        {"gene_index": 1, "name": "tptA", "role": "flanking_conserved", "present": True, "segment_index": 0,
         "start": 120345, "end": 121000, "strand": "+"},
        {"gene_index": 2, "name": "sexM", "role": "core_MAT", "present": False, "segment_index": 0,
         "start": None, "end": None, "strand": None},
    ],
}


def test_write_gff3_includes_only_present_genes(tmp_path):
    out_path = tmp_path / "locus.gff3"
    write_gff3(RECORD, out_path)
    content = out_path.read_text()
    assert "sexP" in content
    assert "tptA" in content
    assert "sexM" not in content
    assert content.startswith("##gff-version 3")


def test_write_proteins_fasta_header_convention(tmp_path):
    out_path = tmp_path / "proteins.faa"
    write_proteins_fasta(RECORD, sequences={0: "MKTAYIAKQRQ", 1: "GATTACAGATTACA"}, out_path=out_path)
    content = out_path.read_text()
    assert ">4837_nrrl-1555_MAT_Plus|gene_index=0|name=sexP|role=core_MAT" in content
    assert "MKTAYIAKQRQ" in content


def test_write_genbank_includes_only_present_genes_and_roundtrips(tmp_path):
    out_path = tmp_path / "locus.gbk"
    write_genbank(RECORD, sequences={0: "MKTAYIAKQRQ", 1: "GATTACAGATTACA"}, out_path=out_path)

    assert out_path.exists()
    content = out_path.read_text()
    assert "sexP" in content
    assert "tptA" in content
    assert "sexM" not in content

    with open(out_path) as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    assert len(records) == 1
    gene_features = {}
    for rec in records:
        for feature in rec.features:
            if feature.type == "gene":
                gene_features[feature.qualifiers["gene"][0]] = feature
    assert set(gene_features) == {"sexP", "tptA"}

    # Pin down the actual coordinate conversion so an off-by-one regression doesn't pass
    # silently. sexP is at absolute 1-based fully-closed [121002, 122400] on a segment
    # starting at absolute 1-based 120345. The placeholder SeqRecord's sequence spans only
    # this segment, so feature locations must be rebased relative to the segment start, not
    # left as absolute coordinates (which would fall outside the record's own sequence bounds
    # and produce invalid GenBank output).
    # relative_start_0based = 121002 - 120345 = 657
    # relative_end_halfopen = 122400 - 120345 + 1 = 2056
    sexp_location = gene_features["sexP"].location
    assert int(sexp_location.start) == 657
    assert int(sexp_location.end) == 2056


# append to tests/db/test_gff_export.py
from unittest.mock import MagicMock


def test_write_genbank_uses_real_sequence_when_ncbi_client_given(tmp_path):
    record = {
        "record_id": "111_a_MAT_combined",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 100, "end": 130,
             "sequence_source": {"type": "insdc_nucleotide", "accession": "ACC1.1", "seq_region": "ACC1.1"}},
        ]}},
        "genes": [
            {"gene_index": 0, "name": "G1", "role": "core_MAT", "present": True,
             "segment_index": 0, "start": 100, "end": 130, "strand": "+"},
        ],
    }
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.return_value = "ATG" * 10 + "TAA"

    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={0: "M" * 10}, out_path=out_path, ncbi=fake_ncbi)

    fake_ncbi.fetch_nucleotide_sequence.assert_called_once_with("ACC1.1", 100, 130, None)
    text = out_path.read_text()
    # Bio.SeqIO's genbank writer always lowercases the ORIGIN sequence block regardless of
    # input case, so sequence-content checks compare case-insensitively.
    assert "N" * 31 not in text.upper()  # the old placeholder is gone
    assert "ATGATGATG" in text.replace("\n", "").replace(" ", "").upper()  # real sequence is present
    assert "CDS" in text
    assert "/translation=" in text.replace("\n", "").replace(" ", "")


def test_write_genbank_falls_back_to_placeholder_when_fetch_fails(tmp_path):
    record = {
        "record_id": "222_b_MAT_combined",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 100, "end": 130,
             "sequence_source": {"type": "insdc_nucleotide", "accession": "ACC2.1", "seq_region": "ACC2.1"}},
        ]}},
        "genes": [],
    }
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.side_effect = Exception("simulated NCBI outage")

    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={}, out_path=out_path, ncbi=fake_ncbi)  # must not raise

    text = out_path.read_text().replace("\n", "").replace(" ", "").upper()
    assert "N" * 31 in text  # placeholder fallback, not a crash


def test_write_genbank_preserves_old_placeholder_behavior_when_no_ncbi_client_given(tmp_path):
    record = {
        "record_id": "333_c_MAT_combined",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 1, "end": 20,
             "sequence_source": {"type": "insdc_nucleotide", "accession": "ACC3.1", "seq_region": "ACC3.1"}},
        ]}},
        "genes": [],
    }
    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={}, out_path=out_path)  # no ncbi= at all -- default behavior

    text = out_path.read_text().replace("\n", "").replace(" ", "").upper()
    assert "N" * 20 in text
