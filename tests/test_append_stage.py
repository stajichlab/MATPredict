"""append_stage.py: a stage's new runs are appended without touching existing rows, and the species that gained runs are decided again with their thin-evidence cap."""
import csv
import importlib.util
import sys
from pathlib import Path

PILOT = Path(__file__).resolve().parents[1] / "datasets/reads/sra_survey/pilot"
spec = importlib.util.spec_from_file_location("append_stage", PILOT / "append_stage.py")
ap = importlib.util.module_from_spec(spec)
sys.modules["append_stage"] = ap
spec.loader.exec_module(ap)

COLS = ["run", "group", "state", "call", "truth", "shared_depth", "screen", "screen_genus", "rescue"]


def _tsv(p, rows, cols):
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(rows)


def _rd(p):
    return list(csv.DictReader(open(p, newline=""), delimiter="\t"))


def test_existing_rows_keep_their_screen_and_rescue_columns():
    final = [{"run": "A", "group": "g", "screen": "mixed", "screen_genus": "Fusarium", "rescue": "genus_panel"}]
    stage = [{"run": "A", "group": "g", "screen": "", "screen_genus": ""}, {"run": "B", "group": "g", "screen": "ok", "screen_genus": "Fusarium"}]
    rows, new = ap.append_rows(final, stage)
    assert rows[0] == final[0] and [r["run"] for r in new] == ["B"] and new[0]["rescue"] == ""


def test_a_topped_up_thin_species_is_decided_again_and_stays_capped(tmp_path):
    w = tmp_path / "wave4"; w.mkdir()
    old = [{"run": f"O{i}", "group": "sp_x", "state": "done", "call": "MAT1-1", "truth": "", "shared_depth": "9.0", "screen": "ok", "screen_genus": "X", "rescue": ""} for i in range(2)]
    other = [{"run": "Z1", "group": "sp_y", "state": "done", "call": "none", "truth": "", "shared_depth": "9.0", "screen": "ok", "screen_genus": "Y", "rescue": "k21"}]
    _tsv(w / "results_all_final.tsv", old + other, COLS)
    gates = [{"slug": "sp_x", "species": "Species x", "decision": "UNVALIDATED:single_locus",
              "reasons": "runs=2 done=2 typed=2 both=0 none/low=0 truth_typed=0 agree=0 median_shared_depth=9.0 thin=single_locus"},
             {"slug": "sp_y", "species": "Species y", "decision": "HOLD:typed_0.00", "reasons": "runs=1 done=1 typed=0 both=0 none/low=1 truth_typed=0 agree=0 median_shared_depth=9.0"}]
    _tsv(w / "gates_final.tsv", gates, ap.GCOLS)
    stage = [{"run": f"N{i}", "group": "sp_x", "state": "done", "call": "MAT1-2", "truth": "MAT1-2", "shared_depth": "9.0", "screen": "ok", "screen_genus": "X"} for i in range(10)]
    _tsv(tmp_path / "stage.tsv", stage, COLS[:-1])
    ap.main(["--wave-dir", str(w), "--stage", str(tmp_path / "stage.tsv")])
    rows = _rd(w / "results_all_final.tsv")
    assert len(rows) == 13 and rows[2]["rescue"] == "k21"                                 # existing rows untouched, 10 appended
    g = {r["slug"]: r for r in _rd(w / "gates_final.tsv")}
    assert g["sp_x"]["decision"] == "UNVALIDATED:single_locus" and "runs=12" in g["sp_x"]["reasons"] and g["sp_x"]["reasons"].endswith("thin=single_locus")
    assert g["sp_y"]["decision"] == "HOLD:typed_0.00"                                       # no new runs: not decided again
    ch = _rd(w / "append_changes.tsv")
    assert [c["slug"] for c in ch] == ["sp_x"]
