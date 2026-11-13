"""Tests for the flank-first scan helpers (scripts/flank_first_scan.py, flank_synteny_score.py)."""
import gzip
import importlib.util
import pathlib
import sys

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ffs = _load("flank_first_scan")
syn = _load("flank_synteny_score")


def test_control_class_ropars_titles():
    assert ffs.control_class("AMN92711.1 hd1-like protein [Rhizophagus irregularis]") == "HD1"
    assert ffs.control_class("AMM63126.1 HD2 AMF mating type protein [Rhizophagus irregularis]") == "HD2"
    assert ffs.control_class("X choline transporter-like protein") == "choline"
    assert ffs.control_class("X phosphoglycerate mutase family protein") == "PGAM"
    assert ffs.control_class("X hypothetical protein") == "other"


def test_control_class_curated_pipe_header():
    assert ffs.control_class("101119_atcc11807_MAT_Minus|gene0|sexM") == "sexM"
    assert ffs.control_class("abc|gene1|rnhA extra") == "rnhA"


def test_parse_gff_orders_genes_and_maps_mrna(tmp_path):
    gff = tmp_path / "a.gff3.gz"
    lines = [
        "c1\tf\tgene\t500\t900\t.\t+\t.\tID=G2;",
        "c1\tf\tmRNA\t500\t900\t.\t+\t.\tID=G2-T1;Parent=G2;",
        "c1\tf\tgene\t10\t100\t.\t-\t.\tID=G1;",
        "c1\tf\tmRNA\t10\t100\t.\t-\t.\tID=G1-T1;Parent=G1;",
    ]
    with gzip.open(gff, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    order, m2g = ffs.parse_gff(str(gff))
    assert [g[0] for g in order] == ["G1", "G2"]
    assert m2g["G2-T1"] == "G2"


def test_read_fasta_plain_and_gz(tmp_path):
    p = tmp_path / "x.faa"
    p.write_text(">a desc\nMK\nLL\n>b\nMM\n")
    assert list(ffs.read_fasta(str(p))) == [("a", "MKLL"), ("b", "MM")]
    z = tmp_path / "x.faa.gz"
    with gzip.open(z, "wt") as fh:
        fh.write(">a\nMK\n")
    assert list(ffs.read_fasta(str(z))) == [("a", "MK")]


def test_block_largest_collinear_set():
    pts = [("c1", 10), ("c1", 12), ("c1", 14), ("c2", 5), ("c1", 40)]
    size, span = syn.block(pts, 5)
    assert size == 3 and span == ("c1", 10, 14)


def test_block_gap_splits_and_empty():
    assert syn.block([("c1", 1), ("c1", 20)], 5)[0] == 1
    assert syn.block([], 5) == (0, None)


def test_n50(tmp_path):
    fai = tmp_path / "g.fai"
    fai.write_text("a\t100\t0\t0\t0\nb\t60\t0\t0\t0\nc\t40\t0\t0\t0\n")
    assert syn.n50(str(fai)) == 100
