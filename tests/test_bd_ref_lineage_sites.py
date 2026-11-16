"""Tests for scripts/bd_ref_lineage_sites.py (lineage-diagnostic divergent sites)."""
import importlib.util
import os

_p = os.path.join(os.path.dirname(__file__), "..", "scripts", "bd_ref_lineage_sites.py")
_s = importlib.util.spec_from_file_location("bd_ref_lineage_sites", _p)
m = importlib.util.module_from_spec(_s)
_s.loader.exec_module(m)


def _paf(tmp_path, name, cs, tstart=0, tend=10, block=10):
    line = "\t".join(["q", "100", "0", "10", "+", "t1", "100", str(tstart), str(tend), "9", str(block), "60", "tp:A:P", f"cs:Z:{cs}"])
    p = tmp_path / name
    p.write_text(line + "\n")
    return str(p)


def test_states(tmp_path):
    # target t1 positions 1-10; substitution at position 4 (A in query)
    paf = _paf(tmp_path, "a.paf", ":3*ga:6", 0, 10)
    subs, blocks = m.parse_paf(paf, min_block=5)
    assert m.site_state(subs, blocks, "t1", 4, "A") == "jel"
    assert m.site_state(subs, blocks, "t1", 4, "C") == "third"
    assert m.site_state(subs, blocks, "t1", 5, "A") == "clft"
    assert m.site_state(subs, blocks, "t1", 50, "A") == "na"


def test_min_block_and_secondary(tmp_path):
    paf = _paf(tmp_path, "b.paf", ":10", 0, 10, block=10)
    subs, blocks = m.parse_paf(paf, min_block=20)
    assert not blocks
    p2 = tmp_path / "c.paf"
    p2.write_text("\t".join(["q", "100", "0", "10", "+", "t1", "100", "0", "10", "10", "10", "60", "tp:A:S", "cs:Z::10"]) + "\n")
    assert not m.parse_paf(str(p2), min_block=5)[1]


def test_two_blocks_cover_is_na(tmp_path):
    p = tmp_path / "d.paf"
    row = "\t".join(["q", "100", "0", "10", "+", "t1", "100", "0", "10", "10", "10", "60", "tp:A:P", "cs:Z::10"])
    p.write_text(row + "\n" + row + "\n")
    subs, blocks = m.parse_paf(str(p), min_block=5)
    assert m.site_state(subs, blocks, "t1", 3, "A") == "na"


def test_is_diagnostic():
    st = {"b1": "clft", "b2": "na", "o1": "jel", "o2": "na"}
    assert m.is_diagnostic(st, ["b1", "b2"], ["o1", "o2"])
    assert not m.is_diagnostic({**st, "b2": "jel"}, ["b1", "b2"], ["o1", "o2"])
    assert not m.is_diagnostic({**st, "o2": "clft"}, ["b1", "b2"], ["o1", "o2"])
    assert not m.is_diagnostic({"b1": "na", "b2": "na", "o1": "jel", "o2": "jel"}, ["b1", "b2"], ["o1", "o2"])
