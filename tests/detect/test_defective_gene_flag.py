"""A polished model with a frameshift or an internal stop codon is flagged as a defective (pseudogenised) gene (issue #151).

Measured case: Lodderomyces parapsilosis has a defective MTLa1 (Logue et al. 2005); with the intact L. orthopsilosis MTLA1 as a
reference, exonerate models the pseudogene and the report called it an ordinary polished gene. The flag is report-only:
the tier and `polished_genes` do not read it.
"""
import dataclasses
import random

from Bio.Seq import Seq

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import (
    DetectionResult, GeneEvidence, LocusSegment, _flag_internal_stops, count_internal_stops,
)
from MATPredict.detect.report import _result_doc, defective_genes

# A long CDS of random sense codons (fixed seed): the true frame has no stop and, at this length, the other two
# frames contain stops, which `count_internal_stops` needs (it takes the frame with the fewest).
_RNG = random.Random(11)
_SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")]
CDS = "ATG" + "".join(_RNG.choice(_SENSE) for _ in range(149)) + "CTGCTG"  # 152 codons incl. two CTG (Ser in table 12)


def _revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _ev(gene="MTLA1", exons=((11, 11 + len(CDS) - 1),), strand="+", status="polished_single", frameshifts=0, **kw):
    return GeneEvidence(gene, "core_MAT", "c1", exons[0][0], exons[-1][1], strand, 48.0, 98.0, "rec", "exonerate_refine",
                        status=status, frameshifts=frameshifts, exons=tuple(exons), **kw)


def _genome(tmp_path, seq):
    p = tmp_path / "g.fa"
    p.write_text(">c1\n" + seq + "\n")
    return p


def test_the_test_cds_reads_through_in_frame_one_and_not_in_the_others():
    for frame in (1, 2):
        prot = str(Seq(CDS[frame:][: (len(CDS) - frame) // 3 * 3]).translate(table=12))
        assert "*" in prot, "the alternate frames must contain a stop for these tests to be meaningful"
    assert "*" not in str(Seq(CDS).translate(table=12))


def test_a_clean_cds_has_no_internal_stop():
    seq = Seq("N" * 10 + CDS + "TAA")
    assert count_internal_stops(seq, [(11, 10 + len(CDS) + 3)], "+", 1) == 0


def test_the_gene_own_trailing_stop_is_not_counted_but_internal_ones_are():
    with_stops = CDS[:30] + "TGA" + CDS[33:60] + "TAG" + CDS[63:]
    assert count_internal_stops(Seq(with_stops), [(1, len(with_stops))], "+", 1) == 2
    assert count_internal_stops(Seq(with_stops + "TAA"), [(1, len(with_stops) + 3)], "+", 1) == 2


def test_the_frame_is_not_assumed_so_a_correct_model_is_never_flagged():
    # the model starts one base into a codon: the true frame is frame 1, which has no stop
    seq = "G" + CDS + "TAA"
    assert count_internal_stops(Seq(seq), [(1, len(seq))], "+", 1) == 0


def test_ctg_is_a_stop_neither_in_the_standard_code_nor_in_table_12():
    assert count_internal_stops(Seq("CTG" * 40), [(1, 120)], "+", 12) == 0
    assert count_internal_stops(Seq("CTG" * 40), [(1, 120)], "+", 1) == 0
    # TGA is a stop in both codes, so it is counted either way
    seq = CDS[:90] + "TGA" + CDS[93:]
    for table in (1, 12):
        assert count_internal_stops(Seq(seq), [(1, len(seq))], "+", table) == 1


def _two_exon_genome(cds, intron="GTAAGTTTTTAG", cut=60):
    """Genome with `cds` split into two exons around `intron`, on the reverse strand: returns (Seq, exon spans)."""
    top = cds[:cut] + intron + cds[cut:]
    genome = "A" * 10 + _revcomp(top) + "A" * 10
    first = len(cds) - cut                     # genomic exon 1 is the 3' CDS part (reverse strand)
    return Seq(genome), [(11, 10 + first), (11 + first + len(intron), 10 + len(top))]


def test_exons_are_spliced_and_the_minus_strand_is_read_from_the_reverse_complement():
    stopped = CDS[:90] + "TAA" + CDS[93:]
    genome, exons = _two_exon_genome(stopped)
    assert count_internal_stops(genome, exons, "-", 1) == 1
    genome_clean, exons_clean = _two_exon_genome(CDS)
    assert count_internal_stops(genome_clean, exons_clean, "-", 1) == 0


def test_flag_marks_a_polished_model_with_an_internal_stop(tmp_path):
    seq = "A" * 10 + CDS[:90] + "TGA" + CDS[93:] + "TAA" + "A" * 10
    fasta = _genome(tmp_path, seq)
    flagged = _flag_internal_stops([_ev(exons=((11, 10 + len(CDS) + 3),))], fasta, 12, {})
    assert flagged[0].internal_stops == 1 and flagged[0].defective


def test_flag_leaves_clean_unpolished_and_frameshifted_entries_alone(tmp_path):
    fasta = _genome(tmp_path, "A" * 10 + CDS + "TAA" + "A" * 10)
    clean = _ev(exons=((11, 10 + len(CDS) + 3),))
    unpolished = GeneEvidence("MTLA1", "core_MAT", "c1", 11, 100, "+", 31.0, 19.0, "rec", "tblastn_genome", status="unpolished")
    shifted = _ev(frameshifts=1)
    out = _flag_internal_stops([clean, unpolished, shifted], fasta, 12, {})
    assert out[0].internal_stops == 0 and not out[0].defective
    assert out[1] is unpolished and not out[1].defective
    # a frameshifted model's exons do not read through, so its stops are not counted; it is defective by the frameshift
    assert out[2].internal_stops == 0 and out[2].defective


def test_an_unreadable_genome_skips_the_flag_instead_of_failing(tmp_path):
    entry = _ev(exons=((11, 10 + len(CDS) + 3),))
    assert _flag_internal_stops([entry], tmp_path / "missing.fa", 12, {}) == [entry]


def test_a_frameshift_alone_marks_the_gene_defective():
    assert _ev(frameshifts=2).defective
    assert not _ev().defective


def _result(genes):
    return DetectionResult(
        family_key=FamilyKey("Ascomycota", "MTL"), contig="c1", start=1, end=5000, confidence="high", idiomorph="A",
        ambiguous_with=[], genes_found=[g.gene_name for g in genes], genes_missing=[], fragmented=False,
        segments=[LocusSegment("c1", 1, 5000, contig_edge_distance=99)], gene_evidence=genes,
        reference_records=["rec"], locus_class="mat_locus", polished_genes=len(genes),
    )


def test_the_report_names_the_defective_genes_and_leaves_tier_and_polished_count_alone():
    intact = _ev("MTLA2")
    pseudo = _ev("MTLA1", frameshifts=1)
    stopped = dataclasses.replace(_ev("PIK1"), internal_stops=2)
    doc = _result_doc(_result([intact, pseudo, stopped]))
    assert doc["defective_genes"] == ["MTLA1", "PIK1"]
    by_gene = {e["gene"]: e for e in doc["gene_evidence"]}
    assert "defective" not in by_gene["MTLA2"] and "internal_stops" not in by_gene["MTLA2"]
    assert by_gene["MTLA1"]["defective"] is True and by_gene["MTLA1"]["frameshifts"] == 1
    assert by_gene["PIK1"]["defective"] is True and by_gene["PIK1"]["internal_stops"] == 2
    assert doc["confidence"] == "high" and doc["polished_genes"] == 3


def test_an_intact_locus_reports_exactly_as_before():
    doc = _result_doc(_result([_ev("MTLA1"), _ev("MTLA2")]))
    assert "defective_genes" not in doc
    assert all("defective" not in e and "internal_stops" not in e for e in doc["gene_evidence"])


def test_the_losing_half_of_a_resolved_pair_is_not_this_loci_defective_gene():
    loser = _ev("MTLA1", frameshifts=1, lost_to="MTLalpha1")
    assert defective_genes(_result([loser, _ev("MTLalpha1")])) == []
