"""Annotation tracks: ORF checks, helixer and GenBank parsing, and the extension summary."""
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq as BioSeq
from Bio.SeqFeature import SeqFeature as F
from Bio import SeqIO

from MATPredict.detect.tracks import (
    Track, TrackGene, coding_sequence, from_genbank, from_helixer_gff3, orf_status, summary,
)

# A 60 bp coding sequence (ATG ... TAA) with no internal stop, on contig "c1".
CDS = "ATGGCTGCTAAAGGTCCAGATTTCTTGGAAGGAGCTGCTAGAAGCCCATATGCAGAAGCTTAA"
CONTIG = "N" * 100 + CDS + "N" * 100


def test_orf_status_complete_orf():
    start, stop, length, complete = orf_status(CDS)
    assert (start, stop, complete) == ("ATG", "TAA", True)
    assert length == len(CDS) // 3


def test_orf_status_flags_missing_start_and_internal_stop():
    assert orf_status("GCTGCTTAA")[3] is False                  # no ATG
    assert orf_status("ATGTAAGCTTAA")[3] is False               # internal stop


def test_helixer_gff_model_exons_and_orf(tmp_path: Path):
    gff = tmp_path / "h.gff3"
    body = CDS
    # split the CDS into two exons, 1-based on c1: 101-120 and 121-160
    last = 100 + len(CDS)
    e1, e2 = (101, 120), (121, last)
    gff.write_text(
        "##gff-version 3\n"
        f"c1\tHelixer\tgene\t101\t{last}\t.\t+\t.\tID=g1\n"
        f"c1\tHelixer\tmRNA\t101\t{last}\t.\t+\t.\tID=g1.t1;Parent=g1\n"
        f"c1\tHelixer\tCDS\t{e1[0]}\t{e1[1]}\t.\t+\t0\tID=cds1;Parent=g1.t1\n"
        f"c1\tHelixer\tCDS\t{e2[0]}\t{e2[1]}\t.\t+\t0\tID=cds2;Parent=g1.t1\n"
    )
    track = from_helixer_gff3(gff, {"c1": CONTIG})
    (g,) = track.genes
    assert g.exons == (e1, e2) and (g.start, g.end, g.strand) == (101, last, "+")
    assert g.source == "helixer" and g.complete_orf is True and g.start_codon == "ATG"


def test_genbank_cds_with_join_and_minus_strand(tmp_path: Path):
    # the same CDS on the minus strand, joined from two parts
    rc = str(Seq(CDS).reverse_complement())
    contig = "N" * 100 + rc + "N" * 100
    rec = SeqRecord(BioSeq(contig), id="c1", name="c1", description="test")
    rec.annotations["molecule_type"] = "DNA"
    loc = CompoundLocation([FeatureLocation(100, 130, strand=-1), FeatureLocation(130, 100 + len(CDS), strand=-1)])
    rec.features.append(SeqFeature(loc, type="CDS", qualifiers={"gene": ["g2"]}))
    path = tmp_path / "a.gbk"
    SeqIO.write(rec, path, "genbank")
    track = from_genbank(path)
    (g,) = track.genes
    assert g.name == "g2" and g.strand == "-" and g.source == "genbank"
    assert g.start == 101 and g.end == 100 + len(CDS)
    assert g.complete_orf is True


def test_summary_extension_signs_on_both_strands():
    ev_plus = Track("evidence", (TrackGene("gp", "c1", 200, 300, "+", ((200, 300),), "evidence"),))
    other_plus = Track("helixer", (TrackGene("hp", "c1", 180, 340, "+", ((180, 340),), "helixer",
                                             complete_orf=True),))
    (row,) = summary(ev_plus, other_plus)
    assert (row.extends_5p, row.extends_3p) == (20, 40)      # starts 20 bp earlier, ends 40 bp later
    assert row.strand_agrees and row.other_complete_orf is True

    ev_minus = Track("evidence", (TrackGene("gm", "c1", 200, 300, "-", ((200, 300),), "evidence"),))
    other_minus = Track("helixer", (TrackGene("hm", "c1", 180, 340, "-", ((180, 340),), "helixer"),))
    (row,) = summary(ev_minus, other_minus)
    # on the minus strand the 5' end is the high coordinate: 340 - 300 = 40; the 3' end is low: 200 - 180 = 20
    assert (row.extends_5p, row.extends_3p) == (40, 20)


def test_summary_reports_other_strand_as_disagreeing_without_extension():
    ev = Track("evidence", (TrackGene("g", "c1", 200, 300, "+", ((200, 300),), "evidence"),))
    other = Track("helixer", (TrackGene("h", "c1", 220, 280, "-", ((220, 280),), "helixer"),))
    (row,) = summary(ev, other)
    assert row.strand_agrees is False and (row.extends_5p, row.extends_3p) == (0, 0)
