"""contamination_table.py: one table of the runs the screen flagged, across committed wave folders, with the species name from the gates."""
import csv
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "datasets/reads/sra_survey/pilot/contamination_table.py"


def _tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(header); w.writerows(rows)


def test_flagged_runs_are_listed_genus_mismatch_first_and_clean_runs_left_out(tmp_path):
    pilot = tmp_path / "pilot"
    for wave in ("wave1", "wave3_x", "notawave"):
        (pilot / wave).mkdir(parents=True)
    hdr = ["run", "group", "strain", "study", "state", "call", "truth", "screen", "screen_genus", "screen_second", "screen_hashes", "screen_explained", "screen_deeper"]
    _tsv(pilot / "wave1" / "results_all_final.tsv", hdr, [
        ["R1", "fusarium_oxysporum", "s1", "P1", "done", "MAT1-1", "", "mixed", "Fusarium", "Meyerozyma", "40", "0.05", ""],
        ["R2", "fusarium_oxysporum", "s2", "P1", "done", "MAT1-2", "", "ok", "Fusarium", "", "", "", ""],
        ["R3", "fusarium_oxysporum", "s3", "P2", "done", "none", "", "genus_mismatch", "Meyerozyma", "", "", "", ""]])
    _tsv(pilot / "wave1" / "gates_final.tsv", ["slug", "species", "decision", "reasons"], [["fusarium_oxysporum", "Fusarium oxysporum", "PASS", ""]])
    _tsv(pilot / "wave3_x" / "results_all_final.tsv", hdr, [["R4", "c_grp", "s4", "P3", "done", "a", "", "borderline", "Cryptococcus", "Lodderomyces", "6", "0.8", "yes"]])
    _tsv(pilot / "wave3_x" / "gates_final.tsv", ["slug", "species", "decision", "reasons"], [["c_grp", "Cryptococcus gattii complex", "UNVALIDATED", ""]])
    _tsv(pilot / "notawave" / "results_all_final.tsv", hdr, [["R5", "x", "", "", "done", "a", "", "mixed", "", "", "", "", ""]])
    out = tmp_path / "flags.tsv"
    subprocess.run([sys.executable, str(SCRIPT), str(pilot), str(out)], check=True, capture_output=True)
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert [r["run"] for r in rows] == ["R3", "R4", "R1"]                     # genus_mismatch first, then by species
    assert rows[0]["species"] == "Fusarium oxysporum" and rows[1]["species"] == "Cryptococcus gattii complex"
    assert rows[2]["second_genus"] == "Meyerozyma" and rows[2]["wave"] == "wave1"
