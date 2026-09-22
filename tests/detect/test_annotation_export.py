"""Turn an annotation GenBank file into the pair `detect` needs.

Promoted from a scratch script on 2026-09-22. It is on the critical path for
every proteome-driven run: the ZygoLife LCG tree has 899 annotated genomes
(844 with a `.proteins.fa`), and BFD carries its own annotations, but NEITHER
writes the coordinates into the protein deflines.

That matters more than it sounds. `search._parse_proteome_location` REQUIRES
`contig:start-end:strand` on every defline and raises `ProteomeDeflineError`
without it -- deliberately, because silently skipping unparseable deflines
would turn a malformed input into a confident "no MAT locus found". The
project handoff records that this exact format mistake cost a whole run.
ZygoLife deflines look like `>EDD06_000001-T1 EDD06_000001`: no coordinates.

The proteome path is worth the trouble: ~19 s/genome against ~4 min for
genome-only, measured on the same code.
"""
from pathlib import Path

import pytest

from MATPredict.detect.annotation_export import convert_genbank


GBK = """LOCUS       scaffold_1               600 bp    DNA     linear   PLN 01-JAN-2026
DEFINITION  Test organism.
ACCESSION   scaffold_1
FEATURES             Location/Qualifiers
     source          1..600
     CDS             complement(30..95)
                     /locus_tag="TEST_000001"
                     /product="first gene"
                     /translation="MKVLAAALWCWSTGVQ"
     CDS             200..262
                     /locus_tag="TEST_000002"
                     /protein_id="XP_000002.1"
                     /translation="MSSNRTFDPQVWKA"
     CDS             400..450
                     /locus_tag="TEST_000003"
ORIGIN
        1 atgaaagttc tggccgccgc cctgtggtgc tggtctaccg gcgtgcagta aatgaaagtt
       61 ctggccgccg ccctgtggtg ctggtctacc ggcgtgcagt aaatgaaagt tctggccgcc
      121 gccctgtggt gctggtctac cggcgtgcag taaatgaaag ttctggccgc cgccctgtgg
      181 tgctggtcta ccggcgtgca gtaaatgaaa gttctggccg ccgccctgtg gtgctggtct
      241 accggcgtgc agtaaatgaa agttctggcc gccgccctgt ggtgctggtc taccggcgtg
      301 cagtaaatga aagttctggc cgccgccctg tggtgctggt ctaccggcgt gcagtaaatg
      361 aaagttctgg ccgccgccct gtggtgctgg tctaccggcg tgcagtaaat gaaagttctg
      421 gccgccgccc tgtggtgctg gtctaccggc gtgcagtaaa tgaaagttct ggccgccgcc
      481 ctgtggtgct ggtctaccgg cgtgcagtaa atgaaagttc tggccgccgc cctgtggtgc
      541 tggtctaccg gcgtgcagta aatgaaagtt ctggccgccg ccctgtggtg ctggtctacc
//
"""


@pytest.fixture
def gbk(tmp_path):
    p = tmp_path / "Test_organism.gbk"
    p.write_text(GBK)
    return p


def _deflines(path):
    return [l[1:].strip() for l in path.read_text().splitlines() if l.startswith(">")]


def test_it_writes_both_files_named_after_the_record(gbk, tmp_path):
    out = convert_genbank(gbk, tmp_path / "out")
    assert out.genome_fasta.name == "Test_organism.fna"
    assert out.proteome_fasta.name == "Test_organism.faa"
    assert out.genome_fasta.exists() and out.proteome_fasta.exists()


def test_every_protein_defline_carries_the_location(gbk, tmp_path):
    """The whole reason this tool exists."""
    from MATPredict.detect.search import _parse_proteome_location
    out = convert_genbank(gbk, tmp_path / "out")
    lines = _deflines(out.proteome_fasta)
    assert lines, "no proteins written"
    for d in lines:
        contig, start, end, strand = _parse_proteome_location(d, d.split()[0])
        assert contig == "scaffold_1"
        assert start < end and strand in "+-"


def test_coordinates_and_strand_are_right(gbk, tmp_path):
    out = convert_genbank(gbk, tmp_path / "out")
    by_id = {d.split()[0]: d for d in _deflines(out.proteome_fasta)}
    assert "scaffold_1:30-95:-" in by_id["TEST_000001"]
    assert "scaffold_1:200-262:+" in by_id["TEST_000002"]


def test_a_cds_with_no_translation_is_skipped(gbk, tmp_path):
    """TEST_000003 has no /translation. Emitting an empty sequence would make
    diamond fail on the whole file."""
    out = convert_genbank(gbk, tmp_path / "out")
    assert "TEST_000003" not in "".join(_deflines(out.proteome_fasta))
    assert out.proteins == 2


def test_the_genome_fasta_keeps_the_contig_name(gbk, tmp_path):
    """Coordinates in the deflines are meaningless if the contig is renamed."""
    out = convert_genbank(gbk, tmp_path / "out")
    names = _deflines(out.genome_fasta)
    assert names and names[0].split()[0] == "scaffold_1"
    assert out.contigs == 1


def test_it_reports_what_it_wrote(gbk, tmp_path):
    out = convert_genbank(gbk, tmp_path / "out")
    assert (out.contigs, out.proteins) == (1, 2)


def test_a_genbank_with_no_cds_raises_rather_than_writing_an_empty_proteome(tmp_path):
    """An empty .faa is the failure mode that looks like success: diamond
    returns no hits and the genome is silently reported as having no locus."""
    p = tmp_path / "Empty.gbk"
    p.write_text(GBK.replace("     CDS", "     misc_feature"))
    with pytest.raises(ValueError, match="no CDS"):
        convert_genbank(p, tmp_path / "out")
