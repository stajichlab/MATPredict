"""merge_screen.py: screen columns join the final table; rescued calls stay; runs the final table lacked are appended."""
import csv
import importlib.util
import sys
from pathlib import Path

PILOT = Path(__file__).resolve().parents[1] / "datasets/reads/sra_survey/pilot"
spec = importlib.util.spec_from_file_location("merge_screen", PILOT / "merge_screen.py")
ms = importlib.util.module_from_spec(spec)
sys.modules["merge_screen"] = ms
spec.loader.exec_module(ms)


def test_screened_runs_gain_the_columns_and_keep_a_rescued_call():
    final = [{"run": "A", "call": "MAT1-1", "rescue": "k21"}, {"run": "B", "call": "none", "rescue": ""}]
    stage = [{"run": "A", "call": "none", "screen": "ok", "screen_genus": "Fusarium"}, {"run": "C", "call": "MAT1-2", "screen": "mixed", "screen_second": "Candidozyma"}]
    rows, added = ms.merge(final, stage)
    a, b, c = rows
    assert (a["call"], a["rescue"], a["screen"], a["screen_genus"]) == ("MAT1-1", "k21", "ok", "Fusarium")     # the stage's untyped k=31 call does not replace the rescue
    assert b["screen"] == "" and b["screen_second"] == ""
    assert added == 1 and (c["run"], c["screen"], c["rescue"]) == ("C", "mixed", "")


def test_the_script_writes_a_table_with_every_screen_column(tmp_path):
    def tsv(p, rows):
        with open(p, "w", newline="") as fh:
            w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    tsv(tmp_path / "f.tsv", [{"run": "A", "call": "MAT1-1", "rescue": ""}])
    tsv(tmp_path / "s.tsv", [{"run": "A", "call": "MAT1-1", "screen": "ok"}])
    ms.main(["--final", str(tmp_path / "f.tsv"), "--stage", str(tmp_path / "s.tsv"), "--out", str(tmp_path / "o.tsv")])
    out = list(csv.DictReader(open(tmp_path / "o.tsv"), delimiter="\t"))
    assert set(ms.SCREEN) <= set(out[0]) and out[0]["screen"] == "ok" and len(out) == 1


def test_a_no_call_run_with_a_species_call_is_ok():
    stage = [{"run": "A", "screen": "no_call"}, {"run": "B", "screen": "no_call"}, {"run": "C", "screen": "mixed"}]
    rows, _ = ms.merge([], stage, species_of=lambda r: "Clonostachys rosea" if r == "A" else "")
    assert [r["screen"] for r in rows] == ["ok", "no_call", "mixed"]
