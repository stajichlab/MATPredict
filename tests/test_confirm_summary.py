"""confirm_summary.py: a flagged set is confirmed only when GmG names the expected second genus in the unmapped reads and the unmapped
fraction is above what the negative controls show."""
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


def _set(out: Path, sid: str, pairs: int, unmapped: int, genus: str | None, primary=None, second=None, aligned="97.5"):
    d = out / sid
    (d / "gmg").mkdir(parents=True)
    (d / "row.tsv").write_text(f"{sid}\t-\t{pairs}\t{unmapped}\t{aligned}\tok\n")
    if genus is not None:
        rec = {"id": sid, "status": "ok", "guess": {"calls": {"genus": {"taxon": genus}, "species": {"taxon": None, "vote": genus + " sp."}}},
               "qc": {"coverage_est": 4.0, "n_hashes": 900, "mixture": {"primary": primary or genus, "suspected": False, "second": {"taxon": second, "hashes": 20} if second else None}}}
        (d / "gmg" / "results.jsonl").write_text(json.dumps(rec) + "\n")


def _run(tmp_path, ids_rows):
    ids = tmp_path / "ids.tsv"
    ids.write_text("".join(f"{a}\t{b}\n" for a, b in ids_rows))
    out_tsv = tmp_path / "table.tsv"
    cs.main(["--out", str(tmp_path / "out"), "--ids", str(ids), "--table", str(out_tsv)])
    return {r["id"]: r for r in csv.DictReader(open(out_tsv), delimiter="\t")}


def test_a_real_contaminant_is_confirmed_and_a_chance_second_genus_is_not(tmp_path):
    out = tmp_path / "out"
    _set(out, "real", 8_000_000, 640_000, "Candidozyma")              # 8% unmapped, and GmG names Candidozyma
    _set(out, "chance", 8_000_000, 96_000, "Gigaspora")                # 1.2% unmapped, named as something else
    _set(out, "neg1", 8_000_000, 120_000, "Escherichia")               # negative controls: 1.5% unmapped
    _set(out, "neg2", 8_000_000, 80_000, "Homo")
    rows = _run(tmp_path, [("real", "Candidozyma"), ("chance", "Candidozyma"), ("neg1", ""), ("neg2", "")])
    assert rows["real"]["verdict"] == "confirmed" and rows["real"]["unmapped_pct"] == "8.0"
    assert rows["chance"]["verdict"] == "not_confirmed"
    assert rows["neg1"]["verdict"] == "control" and rows["neg2"]["verdict"] == "control"


def test_naming_the_genus_is_not_enough_when_the_unmapped_fraction_is_no_higher_than_the_controls(tmp_path):
    out = tmp_path / "out"
    _set(out, "low", 8_000_000, 100_000, "Malassezia")                 # 1.25%: named, but a negative control has 1.5%
    _set(out, "neg", 8_000_000, 120_000, "Escherichia")
    assert _run(tmp_path, [("low", "Malassezia"), ("neg", "")])["low"]["verdict"] == "not_confirmed"


def test_a_mixture_member_counts_as_naming_the_genus(tmp_path):
    out = tmp_path / "out"
    _set(out, "mix", 8_000_000, 900_000, "Lecanicillium", primary="Lecanicillium", second="Histoplasma")
    _set(out, "neg", 8_000_000, 100_000, "Homo")
    assert _run(tmp_path, [("mix", "Histoplasma"), ("neg", "")])["mix"]["verdict"] == "confirmed"


def test_too_few_unmapped_reads_to_classify_and_a_set_that_did_not_run(tmp_path):
    out = tmp_path / "out"
    _set(out, "empty", 8_000_000, 300, None)                            # no GmG record: too few reads left
    _set(out, "neg", 8_000_000, 100_000, "Homo")
    rows = _run(tmp_path, [("empty", "Candidozyma"), ("neg", ""), ("missing", "Aspergillus")])
    assert rows["empty"]["verdict"] == "no_unmapped" and rows["missing"]["verdict"] == "not_run"


def test_without_negative_controls_the_unmapped_fraction_is_not_required(tmp_path):
    out = tmp_path / "out"
    _set(out, "a", 8_000_000, 50_000, "Aspergillus")
    assert _run(tmp_path, [("a", "Aspergillus")])["a"]["verdict"] == "confirmed"
