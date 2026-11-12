"""datasets/reads/sketch_store/strain_starship_table.py: per-strain Starship table with the bafA class; local isolates counted, never named."""
import csv
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("strain_starship_table", ROOT / "datasets" / "reads" / "sketch_store" / "strain_starship_table.py")
t = importlib.util.module_from_spec(spec); spec.loader.exec_module(t)


def write(path, cols, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)


def test_table(tmp_path):
    cols = ["sample_id", "strain", "source", "excluded", "clone", "lineage", "bafA_frac", "Janus-h4_frac", "Lamia-h4_frac", "bafA_call", "Janus-h4_call", "Lamia-h4_call"]
    rows = [dict(sample_id="S1", strain="eAF1", source="public", excluded="", clone="1", lineage="150", bafA_frac="1", **{"Janus-h4_frac": "0.6", "Lamia-h4_frac": "0.9"}, bafA_call="P", **{"Janus-h4_call": "P", "Lamia-h4_call": "P"}),
            dict(sample_id="S2", strain="eAF2", source="public", excluded="", clone="2", lineage="150", bafA_frac="0", **{"Janus-h4_frac": "0.9", "Lamia-h4_frac": "0"}, bafA_call="-", **{"Janus-h4_call": "P", "Lamia-h4_call": "-"}),
            dict(sample_id="S3", strain="eAF3", source="public", excluded="", clone="3", lineage="134", bafA_frac="0", **{"Janus-h4_frac": "0", "Lamia-h4_frac": "0"}, bafA_call="-", **{"Janus-h4_call": "-", "Lamia-h4_call": "-"}),
            dict(sample_id="local__eAFSECRET", strain="eAFSECRET", source="local", excluded="", clone="4", lineage="134", bafA_frac="1", **{"Janus-h4_frac": "0", "Lamia-h4_frac": "0"}, bafA_call="P", **{"Janus-h4_call": "-", "Lamia-h4_call": "-"}),
            dict(sample_id="S9", strain="AF293", source="public", excluded="", clone="9", lineage="1", bafA_frac="1", **{"Janus-h4_frac": "0", "Lamia-h4_frac": "0"}, bafA_call="P", **{"Janus-h4_call": "-", "Lamia-h4_call": "-"})]
    write(tmp_path / "calls.tsv", cols, rows)
    write(tmp_path / "screen.tsv", ["sample_id", "target", "hits"],
          [dict(sample_id=s, target="bafA_all", hits=h) for s, h in (("S1", 3), ("S2", 4), ("S3", 0), ("local__eAFSECRET", 5), ("S9", 3))])
    t.main(["--calls", str(tmp_path / "calls.tsv"), "--baf-screen", str(tmp_path / "screen.tsv"), "--prefix", "eAF", "--set", "Janus-h4=0.7", "--out", str(tmp_path / "o")])
    calls = {r["strain"]: r for r in csv.DictReader(open(tmp_path / "o" / "calls.tsv"), delimiter="\t")}
    assert set(calls) == {"eAF1", "eAF2", "eAF3"}                                   # prefix filter; local isolate and AF293 not listed
    assert calls["eAF1"]["bafA_class"] == "cluster_linked" and calls["eAF2"]["bafA_class"] == "standalone_divergent" and calls["eAF3"]["bafA_class"] == "none"
    assert calls["eAF1"]["Janus-h4"] == "p" and calls["eAF2"]["Janus-h4"] == "P"      # tuned cut 0.7 demotes 0.6 to partial
    assert calls["eAF1"]["present_elements"] == "bafA;Lamia-h4" and calls["eAF1"]["n_present"] == "2"
    summ = {r["element"]: r for r in csv.DictReader(open(tmp_path / "o" / "summary.tsv"), delimiter="\t")}
    assert summ["bafA_all:cluster_linked"]["isolates_present"] == "2" and summ["bafA_all:standalone_divergent"]["isolates_present"] == "1"   # the local isolate is counted
    md = (tmp_path / "o" / "table.md").read_text()
    assert "SECRET" not in md and "local__" not in md and "4 isolates (1 local, counted but not listed)" in md
    assert "SECRET" not in (tmp_path / "o" / "calls.tsv").read_text() and "SECRET" not in (tmp_path / "o" / "fractions.tsv").read_text()
