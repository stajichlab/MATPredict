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
    # The fetched sequence's length must equal the segment's own span (33 nt here);
    # write_genbank rejects a length mismatch as a clamped/failed fetch.
    record = {
        "record_id": "111_a_MAT_combined",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 100, "end": 132,
             "sequence_source": {"type": "insdc_nucleotide", "accession": "ACC1.1", "seq_region": "ACC1.1"}},
        ]}},
        "genes": [
            {"gene_index": 0, "name": "G1", "role": "core_MAT", "present": True,
             "segment_index": 0, "start": 100, "end": 132, "strand": "+"},
        ],
    }
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.return_value = "ATG" * 10 + "TAA"

    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={0: "M" * 10}, out_path=out_path, ncbi=fake_ncbi)

    fake_ncbi.fetch_nucleotide_sequence.assert_called_once_with("ACC1.1", 100, 132, None)
    text = out_path.read_text()
    # Bio.SeqIO's genbank writer always lowercases the ORIGIN sequence block regardless of
    # input case, so sequence-content checks compare case-insensitively.
    assert "N" * 33 not in text.upper()  # the old placeholder is gone
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


# --- assembly-typed segments fetch by seq_region (the coordinates' real reference) ---

def test_write_genbank_assembly_source_fetches_with_seq_region(tmp_path):
    """An `assembly`-typed segment cites a GCA_/GCF_ accession, which efetch cannot
    subrange-fetch. Its `start`/`end` are relative to `seq_region` (the contig), so
    `seq_region` is what gets fetched. Real values from the curated DB:
    `5334_h4-8_Aalpha_4` cites assembly `GCF_000143185.2` with seq_region
    `NW_026089539.1`."""
    record = {
        "record_id": "444_d_MAT_combined",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 100, "end": 130,
             "sequence_source": {"type": "assembly", "accession": "GCF_000143185.2",
                                 "seq_region": "NW_026089539.1"}},
        ]}},
        "genes": [],
    }
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.return_value = "ACGT" * 7 + "TAA"

    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={}, out_path=out_path, ncbi=fake_ncbi)

    fake_ncbi.fetch_nucleotide_sequence.assert_called_once_with("NW_026089539.1", 100, 130, None)
    text = out_path.read_text().replace("\n", "").replace(" ", "").upper()
    assert "N" * 31 not in text
    assert "ACGTACGTACGT" in text


def test_write_genbank_insdc_nucleotide_source_still_fetches_with_accession(tmp_path):
    """Regression guard for the widened branch: an `insdc_nucleotide` segment must
    keep fetching with `accession`, not switch to `seq_region`."""
    record = {
        "record_id": "555_e_MAT_combined",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 10, "end": 40,
             "sequence_source": {"type": "insdc_nucleotide", "accession": "ACC9.1",
                                 "seq_region": "NOT_THIS_ONE.1"}},
        ]}},
        "genes": [],
    }
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.return_value = "GGGG" * 7 + "TAA"

    write_genbank(record, sequences={}, out_path=tmp_path / "locus.gbk", ncbi=fake_ncbi)

    fake_ncbi.fetch_nucleotide_sequence.assert_called_once_with("ACC9.1", 10, 40, None)


def test_write_genbank_assembly_source_falls_back_to_placeholder_when_fetch_fails(tmp_path):
    record = {
        "record_id": "666_f_MAT_combined",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 100, "end": 130,
             "sequence_source": {"type": "assembly", "accession": "GCA_000000000.1",
                                 "seq_region": "JAAGWA010000001.1"}},
        ]}},
        "genes": [],
    }
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.side_effect = Exception("simulated NCBI outage")

    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={}, out_path=out_path, ncbi=fake_ncbi)  # must not raise

    text = out_path.read_text().replace("\n", "").replace(" ", "").upper()
    assert "N" * 31 in text


# --- reading frame / genetic code qualifiers, and the segment-length guard ---
import warnings

import pytest


def _frame_record(gene_extra: dict, start: int, end: int) -> dict:
    gene = {
        "gene_index": 0, "name": "G1", "role": "core_MAT", "present": True,
        "segment_index": 0, "start": start, "end": end, "strand": "+",
    }
    gene.update(gene_extra)
    return {
        "record_id": "999_frame_MAT_test",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": start, "end": end,
             "sequence_source": {"type": "insdc_nucleotide", "accession": "ACC9.1", "seq_region": "ACC9.1"}},
        ]}},
        "genes": [gene],
    }


def _only_cds(gbk_path):
    records = list(SeqIO.parse(gbk_path, "genbank"))
    cds = [f for rec in records for f in rec.features if f.type == "CDS"]
    assert len(cds) == 1
    return records, cds[0]


def test_write_genbank_emits_codon_start_when_gene_is_not_in_frame_one(tmp_path):
    # 2 leading bases, then ATG GCT TAA -> "MA" under codon_start=3.
    nucleotides = "GG" + "ATGGCTTAA"
    record = _frame_record({"codon_start": 3}, 1, len(nucleotides))
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.return_value = nucleotides

    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={0: "MA"}, out_path=out_path, ncbi=fake_ncbi)

    _, cds = _only_cds(out_path)
    assert cds.qualifiers["codon_start"] == ["3"]


def test_write_genbank_emits_transl_table_when_gene_uses_an_alternative_code(tmp_path):
    record = _frame_record({"transl_table": 12}, 1, 9)
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.return_value = "ATGCTGTAA"

    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={0: "MS"}, out_path=out_path, ncbi=fake_ncbi)

    _, cds = _only_cds(out_path)
    assert cds.qualifiers["transl_table"] == ["12"]


def test_write_genbank_omits_codon_start_and_transl_table_at_their_genbank_defaults(tmp_path):
    # Deliberate: GenBank's default for both qualifiers is 1, so a default-valued gene
    # writes neither, keeping the conditional style of gene_class/present_in_idiomorphs.
    record = _frame_record({"codon_start": 1, "transl_table": 1}, 1, 9)
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.return_value = "ATGGCTTAA"

    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={0: "MA"}, out_path=out_path, ncbi=fake_ncbi)

    _, cds = _only_cds(out_path)
    assert "codon_start" not in cds.qualifiers
    assert "transl_table" not in cds.qualifiers


def _translate_written_cds(seq_record, cds) -> str:
    """Translate a WRITTEN CDS feature exactly as a consumer of the file would:
    extract its location from the record's own sequence, honour the written
    /codon_start (1-based, relative to the first base of the CDS) and
    /transl_table (defaulting to 1 when absent, as GenBank defines)."""
    nucleotides = str(cds.extract(seq_record.seq))
    offset = int(cds.qualifiers.get("codon_start", ["1"])[0]) - 1
    table = int(cds.qualifiers.get("transl_table", ["1"])[0])
    framed = nucleotides[offset:]
    framed = framed[: len(framed) - len(framed) % 3]
    return str(Seq(framed).translate(table=table)).rstrip("*")


@pytest.mark.parametrize(
    ("gene_extra", "nucleotides", "protein"),
    [
        # plain, in-frame, standard code
        ({}, "ATGGCTTAA", "MA"),
        # codon_start=3: the first two bases are not part of the reading frame
        ({"codon_start": 3}, "GGATGGCTTAA", "MA"),
        # transl_table=12 (alternative yeast nuclear): CTG is Ser, not Leu
        ({"transl_table": 12}, "ATGCTGTAA", "MS"),
    ],
)
def test_written_cds_round_trips_to_its_own_translation_qualifier(
    tmp_path, gene_extra, nucleotides, protein
):
    """The regression test for the dropped-qualifier bug: translating the CDS as
    WRITTEN must reproduce the CDS's own /translation qualifier."""
    record = _frame_record(gene_extra, 1, len(nucleotides))
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.return_value = nucleotides

    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={0: protein}, out_path=out_path, ncbi=fake_ncbi)

    records, cds = _only_cds(out_path)
    assert _translate_written_cds(records[0], cds) == cds.qualifiers["translation"][0]
    assert cds.qualifiers["translation"][0] == protein


def test_round_trip_helper_fails_when_the_frame_qualifiers_are_ignored(tmp_path):
    """Proof the round-trip test above can actually fail: translating the same
    codon_start=3 CDS in frame 1 does NOT give its /translation."""
    nucleotides = "GGATGGCTTAA"
    record = _frame_record({"codon_start": 3}, 1, len(nucleotides))
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.return_value = nucleotides

    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={0: "MA"}, out_path=out_path, ncbi=fake_ncbi)

    records, cds = _only_cds(out_path)
    naive = str(Seq(str(cds.extract(records[0].seq))[:9]).translate()).rstrip("*")
    assert naive != cds.qualifiers["translation"][0]


def test_write_genbank_rejects_a_fetch_whose_length_differs_from_the_segment(tmp_path):
    """A clamped/truncated efetch must fall back to the all-N placeholder (never pad,
    never fabricate) and warn, instead of writing features past the sequence end."""
    record = {
        "record_id": "333_c_MAT_combined",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 100, "end": 130,
             "sequence_source": {"type": "assembly", "accession": "GCA_1.1", "seq_region": "CONTIG1.1"}},
        ]}},
        "genes": [
            {"gene_index": 0, "name": "G1", "role": "core_MAT", "present": True,
             "segment_index": 0, "start": 100, "end": 130, "strand": "+"},
        ],
    }
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.return_value = "ACGT" * 3  # 12 nt, segment is 31 nt

    out_path = tmp_path / "locus.gbk"
    with pytest.warns(UserWarning, match="12 nt"):
        write_genbank(record, sequences={0: "MMMM"}, out_path=out_path, ncbi=fake_ncbi)

    seq_record = next(SeqIO.parse(out_path, "genbank"))
    assert str(seq_record.seq) == "N" * 31  # placeholder, not the short fetch, not padded
    assert max(int(f.location.end) for f in seq_record.features) <= len(seq_record.seq)


def test_write_genbank_accepts_a_fetch_whose_length_matches_the_segment(tmp_path):
    record = {
        "record_id": "444_d_MAT_combined",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 100, "end": 130,
             "sequence_source": {"type": "insdc_nucleotide", "accession": "ACC4.1", "seq_region": "ACC4.1"}},
        ]}},
        "genes": [],
    }
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.return_value = "A" * 31

    out_path = tmp_path / "locus.gbk"
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # no warning may be raised for a correct length
        write_genbank(record, sequences={}, out_path=out_path, ncbi=fake_ncbi)

    seq_record = next(SeqIO.parse(out_path, "genbank"))
    assert str(seq_record.seq).upper() == "A" * 31
