"""Tests for scripts/tx_six_frame.py."""
import importlib.util
import os

_p = os.path.join(os.path.dirname(__file__), "..", "scripts", "tx_six_frame.py")
_s = importlib.util.spec_from_file_location("tx_six_frame", _p)
m = importlib.util.module_from_spec(_s)
_s.loader.exec_module(m)


def test_translate_and_revcomp():
    assert m.translate("ATGGCCTAA") == "MA*"
    assert m.revcomp("ATGC") == "GCAT"
    assert m.translate("ATGNNN") == "MX"


def test_main_writes_six_frames(tmp_path, monkeypatch):
    fa = tmp_path / "t.fa"
    fa.write_text(">t1 x\nATGGCCTAAGGG\n")
    out = tmp_path / "o.faa"
    monkeypatch.setattr("sys.argv", ["x", "--fasta", str(fa), "--out", str(out)])
    m.main()
    recs = [l for l in out.read_text().splitlines() if l.startswith(">")]
    assert recs == [f">t1|{d}|f{i}" for d in ("sense", "anti") for i in (1, 2, 3)]
    assert out.read_text().splitlines()[1] == "MA*G"
