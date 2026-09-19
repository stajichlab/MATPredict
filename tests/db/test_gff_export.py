from __future__ import annotations

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation

from MATPredict.db.gff_export import _cds_location, write_genbank, write_gff3, write_proteins_fasta

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


# --- CompoundLocation exon-structure tests (real bug fix: write_genbank must build a
# real multi-part CDS location from gene["exons"], not a single genomic-span-including-
# introns location; see _cds_location's docstring for the exon-order convention this
# relies on) ---

_SEGMENT = {"segment_index": 0, "start": 1000, "end": 2000}


def test_cds_location_plus_strand_multi_exon_builds_compound_location_in_stored_order():
    # Plus-strand curated-record convention: exons stored in ascending genomic order,
    # which is already transcript (5'->3') order for a plus-strand gene.
    gene = {
        "start": 1010, "end": 1200, "strand": "+",
        "exons": [{"start": 1010, "end": 1050}, {"start": 1100, "end": 1150}, {"start": 1180, "end": 1200}],
    }
    location = _cds_location(gene, _SEGMENT, strand=1)
    assert isinstance(location, CompoundLocation)
    # relative_start = exon.start - segment.start ; relative_end = exon.end - segment.start + 1
    starts = [int(part.start) for part in location.parts]
    ends = [int(part.end) for part in location.parts]
    assert starts == [10, 100, 180]
    assert ends == [51, 151, 201]
    for part in location.parts:
        assert part.strand == 1


def test_cds_location_minus_strand_compound_location_part_order():
    """The real bug this task's brief warned about: for a minus-strand multi-exon gene,
    the CompoundLocation's `parts` must be in the order that reconstructs the correct
    5'->3' transcript when Biopython's own `.extract()`/GenBank writer walks them --
    NOT simply the flat genomic span, and not naively reversed either. This project's
    curated-record schema already stores `gene["exons"]` in descending genomic
    (transcript) order for a minus-strand gene (verified against the real COX13 record
    in db/Ascomycota/Onygenales/199306_rmscc1040_MAT_MAT1-1/metadata.yaml), so
    `_cds_location` must NOT reverse it again -- reversing would produce ascending
    (wrong) order. This test proves the resulting CompoundLocation actually splices to
    the correct real transcript, not just that the part list order matches input order.
    """
    # Two distinct 10nt blocks on a 40nt segment-relative sequence so a wrong splice
    # order/direction is unambiguous: pos[0:10]="A"*10, pos[20:30]="T"*10.
    # Real transcript (5'->3', minus strand) must read the HIGH-coordinate exon (T block,
    # genomic descending) first, each individually reverse-complemented: revcomp(T*10)
    # + revcomp(A*10) = "A"*10 + "T"*10.
    segment_seq = Seq("A" * 10 + "N" * 10 + "T" * 10 + "N" * 10)
    segment = {"segment_index": 0, "start": 1, "end": 40}
    # Stored in descending genomic order, matching this project's real minus-strand
    # curated-record convention (COX13-shaped): high exon first, low exon second.
    gene = {
        "start": 1, "end": 30, "strand": "-",
        "exons": [{"start": 21, "end": 30}, {"start": 1, "end": 10}],
    }
    location = _cds_location(gene, segment, strand=-1)
    assert isinstance(location, CompoundLocation)
    assert [int(p.start) for p in location.parts] == [20, 0]
    assert [int(p.end) for p in location.parts] == [30, 10]
    for part in location.parts:
        assert part.strand == -1

    spliced = str(location.extract(segment_seq))
    assert spliced == "A" * 10 + "T" * 10


def test_cds_location_falls_back_to_single_span_when_no_exons_recorded():
    gene = {"start": 1010, "end": 1200, "strand": "+", "exons": []}
    location = _cds_location(gene, _SEGMENT, strand=1)
    assert not isinstance(location, CompoundLocation)
    assert int(location.start) == 10
    assert int(location.end) == 201


def test_write_genbank_builds_compound_location_for_real_multi_exon_minus_strand_gene(tmp_path):
    """End-to-end: write_genbank must emit a real CompoundLocation CDS (visible as
    complement(join(...)) in the GenBank text) for a curated multi-exon minus-strand
    gene, not one unbroken span covering introns."""
    record = {
        "record_id": "444_d_MAT_combined",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 1, "end": 40, "sequence_source": {"seq_region": "scaffold_x"}},
        ]}},
        "genes": [
            {"gene_index": 0, "name": "cox13like", "role": "flanking_conserved", "present": True,
             "segment_index": 0, "start": 1, "end": 30, "strand": "-",
             "exons": [{"start": 21, "end": 30}, {"start": 1, "end": 10}]},
        ],
    }
    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={0: "M" * 6}, out_path=out_path)

    text = out_path.read_text()
    assert "join(" in text
    assert "complement(" in text

    with open(out_path) as handle:
        records = list(SeqIO.parse(handle, "genbank"))
    cds_features = [f for f in records[0].features if f.type == "CDS"]
    assert len(cds_features) == 1
    assert isinstance(cds_features[0].location, CompoundLocation)
    assert len(cds_features[0].location.parts) == 2
