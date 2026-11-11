"""datasets/reads/sketch_store/starship_report_data.py: screen summary parsing, heatmap order, 2x2 table, JSON and HTML output."""
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
T = ROOT / "tests"


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


rd = load("starship_report_data", "datasets/reads/sketch_store/starship_report_data.py")
sd = load("starship_dynamics_for_report", "datasets/reads/sketch_store/starship_dynamics.py")
fx = load("test_sd_fixture", "tests/test_sketch_store_starship_dynamics.py")

BEFORE = "bafA                                     len     638  specific hashes      3  present 407  partial 782  absent 3660  | vs truth TP 17 FN 2 FP 9 TN 198\n"
AFTER = ("bafA                                     len     638  specific hashes      2  present 411  partial 0  absent 4438  | vs truth TP 17 FN 2 FP 9 TN 198"
         "  | dropped 1 hash(es) seen alone in 778 samples\nA  len 10  specific hashes 4  present 1  partial 2  absent 3\n")


def test_screen_summary(tmp_path):
    (tmp_path / "b.txt").write_text(BEFORE); (tmp_path / "a.txt").write_text(AFTER)
    b, a = rd.screen_summary(tmp_path / "b.txt"), rd.screen_summary(tmp_path / "a.txt")
    assert b["bafA"]["partial"] == 782 and b["bafA"]["truth"] == dict(TP=17, FN=2, FP=9, TN=198) and b["bafA"]["dropped"] == ""
    assert a["bafA"]["partial"] == 0 and "778" in a["bafA"]["dropped"]
    assert a["A"]["truth"] is None and a["A"]["absent"] == 3


def test_report_data_and_html(tmp_path):
    store = fx.make_store(tmp_path)
    out = tmp_path / "dyn"
    sd.main(["--store", str(store), "--out", str(out), "--min-group", "1", "--truth", str(tmp_path / "truth.tsv")])
    (tmp_path / "b.txt").write_text(BEFORE); (tmp_path / "a.txt").write_text(AFTER)
    (tmp_path / "cand.tsv").write_text("tier\tsample_id\tstrain\nstrong\ts1\tst1\n")
    (tmp_path / "t.html").write_text("<script>const D=/*__DATA__*/null;</script>")
    rd.main(["--default", str(out), "--tuned", str(out), "--candidates", str(tmp_path / "cand.tsv"), "--before", str(tmp_path / "b.txt"),
             "--after", str(tmp_path / "a.txt"), "--pair", "A,B", "--meta", "commit=abc123", "--out", str(tmp_path / "d.json"),
             "--template", str(tmp_path / "t.html"), "--html", str(tmp_path / "r.html")])
    d = json.loads((tmp_path / "d.json").read_text())
    assert d["elements"] == ["A", "B"] and d["meta"]["commit"] == "abc123"
    assert [h["clone"] for h in d["heat"]][:2] == ["c1", "c2"]                    # tree order: t1 (c1), t3 (c2), ...
    assert d["heat"][0]["bits"] == "11" and d["heat"][-1]["clone"] == "c5" and d["heat"][-1]["bits"] == "20"   # 2 = missing
    assert d["baf_neb"] == {"11": 2, "10": 0, "01": 0, "00": 2}                 # c1, c2 have both; c3, c4 neither; c5 is missing A
    assert d["sweep"]["A"][4]["present_cut"] == 0.5 and d["sweep"]["A"][4]["TP"] == 3
    assert d["candidates"][0]["sample_id"] == "s1" and d["rescreen"]["after"]["bafA"]["present"] == 411
    html = (tmp_path / "r.html").read_text()
    assert "/*__DATA__*/" not in html and '"elements":["A","B"]' in html
