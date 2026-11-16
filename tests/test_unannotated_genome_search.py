"""Tests for scripts/unannotated_genome_search.py (query classes, region merging, pairing)."""
import importlib.util
import pathlib

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
spec = importlib.util.spec_from_file_location("unannotated_genome_search", SCRIPTS / "unannotated_genome_search.py")
ug = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ug)


def test_query_class_header_styles():
    assert ug.query_class("3231046_uomd23-19_MAT_Plus|gene1|sexP") == "HMG"
    assert ug.query_class("4837_nrrl1555_MAT_Minus|gene1|sexM") == "HMG"
    assert ug.query_class("5011_unknown-1_PM_combined|gene_index=1|name=matMc|role=core_MAT") == "HMG"
    assert ug.query_class("5011_unknown-1_PM_combined|gene_index=0|name=matPi|role=core_MAT") == "HD"
    assert ug.query_class("4896_972h_mat1_M|gene_index=2|name=mat1-Pc|role=core_MAT") == "HD"
    assert ug.query_class("AMN92711.1 hd1-like protein [Rhizophagus irregularis]") == "HD"
    assert ug.query_class("AMM63126.1 HD2 AMF mating type protein") == "HD"
    assert ug.query_class("3108396_ai-2023a_MAT_Plus|gene1|tptA") is None
    assert ug.query_class("AMM63120.1 choline transporter-like protein") is None


def test_merge_regions():
    hits = [("c1", 100, 300, 50.0, 1e-9, 60.0), ("c1", 1500, 1800, 60.0, 1e-20, 90.0), ("c1", 9000, 9100, 40.0, 1e-6, 30.0),
            ("c2", 500, 400, 45.0, 1e-8, 70.0)]
    reg = ug.merge(hits, 2000)
    assert [(r[0], r[1], r[2]) for r in reg] == [("c1", 100, 1800), ("c1", 9000, 9100), ("c2", 400, 500)]
    assert reg[0][3] == 60.0 and reg[0][5] == 150.0


def test_pairs_distance_and_contig():
    hmg = [("c1", 100, 300, 50.0, 1e-9, 60.0), ("c2", 100, 300, 50.0, 1e-9, 60.0)]
    hd = [("c1", 15000, 15500, 40.0, 1e-9, 50.0), ("c1", 90000, 90100, 40.0, 1e-9, 50.0), ("c3", 100, 200, 40.0, 1e-9, 50.0)]
    assert ug.pairs(hmg, hd, 20000) == [(hmg[0], hd[0])]
    assert ug.pairs(hmg, hd, 100) == []
