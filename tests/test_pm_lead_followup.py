"""Tests for scripts/pm_lead_followup.py (gene lookup and layout helpers)."""
import gzip
import importlib.util
import pathlib

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
spec = importlib.util.spec_from_file_location("pm_lead_followup", SCRIPTS / "pm_lead_followup.py")
pm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pm)

GENES = {"c1": [(100, 200, "+", "g1"), (300, 400, "-", "g2"), (500, 600, "+", "g3")]}


def test_overlapping_gene_picks_largest_overlap():
    assert pm.overlapping_gene(GENES, "c1", 350, 520)[1][3] == "g2"  # 50 bp on g2, 21 bp on g3
    assert pm.overlapping_gene(GENES, "c1", 510, 700)[1][3] == "g3"


def test_overlapping_gene_none_when_no_overlap_or_contig():
    assert pm.overlapping_gene(GENES, "c1", 210, 290) == (None, None)
    assert pm.overlapping_gene(GENES, "c9", 100, 200) == (None, None)


def test_read_gff_genes_sorted_and_gene_rows_only(tmp_path):
    p = tmp_path / "x.gff3.gz"
    with gzip.open(p, "wt") as fo:
        fo.write("##gff-version 3\n")
        fo.write("c1\tf\tgene\t300\t400\t.\t-\t.\tID=g2;Name=b\n")
        fo.write("c1\tf\tmRNA\t300\t400\t.\t-\t.\tID=g2-T1\n")
        fo.write("c1\tf\tgene\t100\t200\t.\t+\t.\tID=g1\n")
    assert pm.read_gff_genes(str(p))["c1"] == [(100, 200, "+", "g1"), (300, 400, "-", "g2")]


def test_layout_row_flank_and_domains():
    hits = {"g2": [("HMG_box", 1e-9), ("Other", 1e-9)]}
    rows = pm.layout_row([1], GENES["c1"], hits, 1)
    assert [r[1] for r in rows] == ["g1", "g2", "g3"]
    assert rows[1][5] == "HMG_box"
    assert pm.layout_row([0], GENES["c1"], {}, 5)[0][0] == 0  # clipped at contig start


def test_augustus_proteins_blocks(tmp_path):
    spec2 = importlib.util.spec_from_file_location("augustus_proteins", SCRIPTS / "augustus_proteins.py")
    ap = importlib.util.module_from_spec(spec2)
    spec2.loader.exec_module(ap)
    lines = ["# start gene g1\n", "c1\tAUGUSTUS\tgene\t1\t90\t.\t+\t.\tID=g1\n",
             "# protein sequence = [MKVL\n", "# AAGG]\n", "# end gene g1\n",
             "# start gene g2\n", "# protein sequence = [MSS]\n", "# end gene g2\n"]
    assert list(ap.proteins(lines)) == [("g1", "MKVLAAGG"), ("g2", "MSS")]
