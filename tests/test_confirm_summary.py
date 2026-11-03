"""confirm_summary.py: a flagged set is confirmed when GmG, run on the reads that do not map to Coccidioides, names the expected genus (as its call, or as
the second organism with its mixture flag raised); a chance second genus without the flag is not enough, and the unmapped fraction does not decide."""
import csv
import importlib.util
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("confirm_summary", REPO / "datasets/reads/contamination_spikein/confirm_summary.py")
cs = importlib.util.module_from_spec(spec)
sys.modules["confirm_summary"] = cs
spec.loader.exec_module(cs)


def _set(out: Path, sid: str, pairs: int, unmapped: int, genus, species=None, primary=None, second=None, mixed=False):
    d = out / sid
    (d / "gmg").mkdir(parents=True)
    (d / "row.tsv").write_text(f"{sid}\t-\t{pairs}\t{unmapped}\t97.5\tok\n")
    if genus is not None or species is not None:
        rec = {"id": sid, "status": "ok", "guess": {"calls": {"genus": {"taxon": genus}, "species": {"taxon": species, "vote": species}}},
               "qc": {"coverage_est": 4.0, "n_hashes": 900, "mixture": {"primary": primary or genus, "suspected": mixed, "second": {"taxon": second, "hashes": 100} if second else None}}}
        (d / "gmg" / "results.jsonl").write_text(json.dumps(rec) + "\n")


def _run(tmp_path, ids_rows, capsys=None):
    ids = tmp_path / "ids.tsv"
    ids.write_text("".join(f"{a}\t{b}\n" for a, b in ids_rows))
    out_tsv = tmp_path / "table.tsv"
    cs.main(["--out", str(tmp_path / "out"), "--ids", str(ids), "--table", str(out_tsv)])
    return {r["id"]: r for r in csv.DictReader(open(out_tsv), delimiter="\t")}


def test_a_flagged_second_genus_in_the_unmapped_reads_is_confirmed_and_a_chance_one_is_not(tmp_path):
    out = tmp_path / "out"
    _set(out, "real", 8_000_000, 640_000, "Coccidioides", second="Candidozyma", mixed=True)       # GmG raises its flag on the unmapped reads
    _set(out, "chance", 8_000_000, 96_000, "Coccidioides", second="Candidozyma", mixed=False)     # named, but no flag: 7 hashes of chance
    _set(out, "neg", 8_000_000, 120_000, "Coccidioides", second="Aspergillus")
    rows = _run(tmp_path, [("real", "Candidozyma"), ("chance", "Candidozyma"), ("neg", "")])
    assert rows["real"]["verdict"] == "confirmed" and rows["chance"]["verdict"] == "not_confirmed" and rows["neg"]["verdict"] == "control"


def test_the_unmapped_fraction_does_not_decide(tmp_path):
    out = tmp_path / "out"
    _set(out, "low", 8_000_000, 100_000, "Coccidioides", second="Malassezia", mixed=True)         # 1.25% unmapped, below a control's 19%: still confirmed
    _set(out, "odd_control", 8_000_000, 1_500_000, "Coccidioides", second="Aspergillus")           # a clean control that maps badly
    rows = _run(tmp_path, [("low", "Malassezia"), ("odd_control", "")])
    assert rows["low"]["verdict"] == "confirmed" and rows["low"]["unmapped_pct"] == "1.25"


def test_the_contaminant_as_the_genus_or_species_call_of_the_unmapped_reads_is_confirmed(tmp_path):
    out = tmp_path / "out"
    _set(out, "hist", 8_000_000, 3_300_000, "Histoplasma", species="Histoplasma capsulatum")      # the unmapped reads are the contaminant
    _set(out, "calonectria", 8_000_000, 480_000, None, species="Calonectria ilicicola", second="Coccidioides")   # no genus call, species named
    _set(out, "wrong", 8_000_000, 480_000, "Fusarium", species="Fusarium oxysporum")
    rows = _run(tmp_path, [("hist", "Histoplasma"), ("calonectria", "Calonectria"), ("wrong", "Histoplasma")])
    assert rows["hist"]["verdict"] == "confirmed" and rows["calonectria"]["verdict"] == "confirmed" and rows["wrong"]["verdict"] == "not_confirmed"


def test_too_few_unmapped_reads_to_classify_and_a_set_that_did_not_run(tmp_path):
    out = tmp_path / "out"
    _set(out, "empty", 8_000_000, 300, None)
    rows = _run(tmp_path, [("empty", "Candidozyma"), ("missing", "Aspergillus")])
    assert rows["empty"]["verdict"] == "no_unmapped" and rows["missing"]["verdict"] == "not_run"


def test_the_summary_reports_how_many_controls_gmg_flagged(tmp_path, capsys):
    out = tmp_path / "out"
    _set(out, "a", 8_000_000, 500_000, "Coccidioides", second="Aspergillus", mixed=True)
    _set(out, "neg1", 8_000_000, 100_000, "Coccidioides", second="Aspergillus", mixed=False)
    _set(out, "neg2", 8_000_000, 100_000, "Coccidioides", second="Malassezia", mixed=True)       # a control GmG flags: worth knowing
    _run(tmp_path, [("a", "Aspergillus"), ("neg1", ""), ("neg2", "")])
    err = capsys.readouterr().err
    assert "1 of 1 flagged sets confirmed" in err and "GmG flagged 1 of 2 negative controls" in err
