"""scripts/zygo_regression.py maps contig calls back to scaffold coordinates."""
import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "zygo_regression", Path(__file__).resolve().parents[1] / "scripts" / "zygo_regression.py")
zr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(zr)


def test_a_forward_component_is_offset_onto_its_scaffold():
    amap = {"contig_552": ("scaffold_145", 31020, 41300, 1, "+")}
    assert zr.to_scaffold("contig_552", 272, 10473, amap) == ("scaffold_145", 31291, 41492)


def test_a_reverse_component_is_flipped():
    amap = {"c9": ("s1", 1001, 2000, 1, "-")}
    assert zr.to_scaffold("c9", 1, 100, amap) == ("s1", 1901, 2000)


def test_an_unmapped_name_is_taken_as_a_scaffold():
    assert zr.to_scaffold("scaffold_7", 5, 9, {}) == ("scaffold_7", 5, 9)


def test_agp_rows_are_read(tmp_path):
    (tmp_path / "x.agp").write_text(
        "scaffold_1\t1\t500\t1\tW\tcontig_1\t1\t500\t+\n"
        "scaffold_1\t501\t600\t2\tN\t100\tscaffold\tyes\tunspecified\n")
    assert zr.agp_map(tmp_path / "x.contigs.fsa") == {"contig_1": ("scaffold_1", 1, 500, 1, "+")}
