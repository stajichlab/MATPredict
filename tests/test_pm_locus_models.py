"""Tests for scripts/pm_locus_models.py (exonerate GFF parsing and CDS translation)."""
import importlib.util
import pathlib

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
import sys
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("pm_locus_models", SCRIPTS / "pm_locus_models.py")
pl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pl)

GFF = """Command line: [exonerate ...]
# --- START OF GFF DUMP ---
#
w\texonerate:protein2genome:local\tgene\t11\t40\t394\t+\t.\tgene_id 0 ; sequence q1 ; gene_orientation +
w\texonerate:protein2genome:local\tcds\t11\t19\t.\t+\t.\t
w\texonerate:protein2genome:local\tcds\t31\t40\t.\t+\t.\t
# --- END OF GFF DUMP ---
# --- START OF GFF DUMP ---
w\texonerate:protein2genome:local\tgene\t1\t9\t100\t-\t.\tgene_id 0 ; sequence q2 ; gene_orientation -
w\texonerate:protein2genome:local\tcds\t1\t9\t.\t-\t.\t
# --- END OF GFF DUMP ---
"""


def test_parse_exonerate_gff_models():
    m = pl.parse_exonerate_gff(GFF)
    assert len(m) == 2
    assert m[0]["score"] == 394.0 and m[0]["query"] == "q1" and m[0]["strand"] == "+"
    assert m[0]["cds"] == [(11, 19), (31, 40)]
    assert m[1]["strand"] == "-" and m[1]["cds"] == [(1, 9)]


def test_cds_protein_spliced_plus_strand():
    # exon 1 = ATGGCCAAA (window bases 11-19), an intron, exon 2 = TAA (bases 31-33): M A K stop
    window = "N" * 10 + "ATGGCCAAA" + "GT" + "A" * 9 + "TAA" + "GGGCCCT"
    assert window[30:33] == "TAA"
    assert pl.cds_protein(window, [(11, 19), (31, 33)], "+") == "MAK"  # the trailing stop is removed


def test_cds_protein_minus_strand():
    from region_pfam_check import revcomp
    dna = "ATGGCCAAA"
    window = "CC" + revcomp(dna) + "GG"
    assert pl.cds_protein(window, [(3, 11)], "-") == "MAK"
