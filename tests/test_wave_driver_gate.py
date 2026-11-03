"""The wave gate: a species with a usable panel is shown flagged (UNVALIDATED) and still gets the full run; only a hard failure holds it."""
import csv
import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PILOT = REPO / "datasets/reads/sra_survey/pilot"
spec = importlib.util.spec_from_file_location("wave_driver", PILOT / "wave_driver.py")
wd = importlib.util.module_from_spec(spec)
sys.modules["wave_driver"] = wd
spec.loader.exec_module(wd)


def decide(done=20, typed=20, both=0, truth=0, agree=0, msd=10.0, counted=None):
    return wd.decide(done if counted is None else counted, done, typed, both, truth, agree, msd)


def word(res):
    return res[0].split(":")[0]


def test_the_strict_cases_are_unchanged():
    assert decide(truth=10, agree=10)[0] == "PASS"                               # all typed, 10 truth runs agree
    assert decide(truth=3, agree=3)[0] == "PROVISIONAL"                          # too few comparisons to confirm
    assert decide(typed=2, both=18)[0] == "PROVISIONAL_HOMOTHALLIC"              # both is the expected answer


def test_typing_between_half_and_the_bar_is_unvalidated_and_below_half_is_held():
    r = decide(typed=14, both=0, truth=10, agree=10)                              # 70% typed
    assert r[0] == "UNVALIDATED:typed_0.70" and r[1] == ["typed_0.70"]
    assert word(decide(typed=10, truth=10, agree=10)) == "UNVALIDATED"            # exactly half still runs in full
    assert word(decide(typed=9, truth=10, agree=10)) == "HOLD"                    # under half does not
    assert decide(typed=15, truth=10, agree=10)[0] == "PASS"                      # 75%: no reason at all


def test_a_provisional_species_under_the_provisional_bar_is_unvalidated_not_held():
    assert decide(typed=16, truth=3, agree=3)[0] == "UNVALIDATED:provisional_typed_0.80"      # Blastomyces-like: 0.81 against 0.85
    assert word(decide(typed=8, truth=3, agree=3)) == "HOLD"


def test_too_many_both_and_low_depth_are_soft_reasons():
    r = decide(typed=14, both=6, truth=10, agree=10)
    assert word(r) == "UNVALIDATED" and "both_0.30" in r[1] and "typed_0.70" in r[1]
    assert decide(msd=2.4, truth=10, agree=10)[0] == "UNVALIDATED:shared_depth_2.4"            # Epichloe-like: only the coverage is short


def test_agreement_holds_only_with_ten_or_more_comparison_runs():
    assert decide(truth=12, agree=10)[0].startswith("HOLD:agree_10of12")           # two disagreements in twelve
    assert decide(truth=7, agree=5)[0] == "UNVALIDATED:agree_5of7"                  # the same shortfall on a small sample is a flag
    assert decide(truth=9, agree=8)[0] == "UNVALIDATED:agree_8of9"                  # 89% is under the 95% bar even with one disagreement: a flag, below ten comparison runs
    assert decide(truth=10, agree=10)[0] == "PASS" and decide(truth=20, agree=19)[0] == "PASS"   # 95% with at most one disagreement passes


def test_failed_downloads_hold_the_species():
    r = decide(done=12, typed=12, counted=20, truth=10, agree=10)                    # 60% finished
    assert r[0].startswith("HOLD:many_failed")
    assert word(decide(done=0, typed=0, counted=0)) == "HOLD"


def test_hard_reasons_come_first_and_soft_ones_stay_in_the_holds_text():
    r = decide(typed=8, both=8, msd=2.0, truth=3, agree=3)
    assert word(r) == "HOLD" and r[1][0].startswith("typed_") and any(x.startswith("shared_depth") for x in r[1])


def _gates(tmp_path, rows):
    d = tmp_path / "wave"
    d.mkdir()
    with open(d / "gates_final.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["slug", "species", "decision", "reasons"])
        w.writerows(rows)
    return d


def _run_regate(d, **kw):
    import argparse
    wd.regate(argparse.Namespace(wave=str(d), in_file="gates_final.tsv", out_file=kw.get("out_file")))


ROWS = [
    ["a", "Species a", "PASS", "runs=50 done=50 typed=49 both=1 none/low=0 truth_typed=25 agree=25 median_shared_depth=13.3"],
    ["b", "Species b", "HOLD:typed_0.70+provisional_typed_0.70", "runs=20 done=20 typed=14 both=0 none/low=6 truth_typed=3 agree=3 median_shared_depth=9.0"],
    ["c", "Species c", "HOLD:typed_0.40", "runs=20 done=20 typed=8 both=0 none/low=12 truth_typed=3 agree=3 median_shared_depth=9.0"],
    ["d", "Species d", "HOLD:panel_build_failed", "held:panel_build_failed"],
    ["e", "Species e", "HOLD:many_failed+typed_0.70", "runs=20 done=12 typed=8 both=0 none/low=4 truth_typed=0 agree=0 median_shared_depth=9.0"],
    ["f", "Species f", "PROVISIONAL", "runs=20 done=20 typed=20 both=0 none/low=0 truth_typed=3 agree=3 median_shared_depth=9.0"],
]


def test_regate_changes_only_the_species_whose_decision_changes_and_records_them(tmp_path):
    d = _gates(tmp_path, ROWS)
    before = (d / "gates_final.tsv").read_text().splitlines()
    _run_regate(d)
    after = (d / "gates_final.tsv").read_text().splitlines()
    assert [i for i, (x, y) in enumerate(zip(before, after)) if x != y] == [2]       # only species b moved (b is line 3: header is line 1)
    assert after[2].split("\t")[2] == "UNVALIDATED:typed_0.70+provisional_typed_0.70"
    changes = list(csv.DictReader(open(d / "regate_changes.tsv"), delimiter="\t"))
    assert [(c["slug"], c["previous"].split(":")[0], c["decision"].split(":")[0]) for c in changes] == [("b", "HOLD", "UNVALIDATED")]
    assert after[3].split("\t")[2] == "HOLD:typed_0.40" and after[4].split("\t")[2] == "HOLD:panel_build_failed"      # held for hard reasons or not tested: untouched
    assert after[5].split("\t")[2].startswith("HOLD:many_failed")                                                       # failed downloads stay held


def test_regate_twice_changes_nothing_the_second_time(tmp_path):
    d = _gates(tmp_path, ROWS)
    _run_regate(d)
    once = (d / "gates_final.tsv").read_text()
    _run_regate(d)
    assert (d / "gates_final.tsv").read_text() == once
    assert list(csv.DictReader(open(d / "regate_changes.tsv"), delimiter="\t")) == []


def test_the_committed_gate_tables_are_what_regate_gives():
    """The committed tables of both waves are already regated: nothing would change, no species is held for a soft reason alone, and
    every species that is not held keeps its decision."""
    for wave in ("wave1", "wave2"):
        for r in csv.DictReader(open(PILOT / wave / "gates_final.tsv", newline=""), delimiter="\t"):
            m = wd.REGATE.search(r["reasons"])
            if not m:
                continue
            done, typed, both, truth, agree, msd = int(m[2]), int(m[3]), int(m[4]), int(m[6]), int(m[7]), float(m[8])
            counted = done if "many_failed" not in r["decision"] else done * 2 + 1
            new, _ = wd.decide(counted, done, typed, both, truth, agree, msd)
            assert new.split(":")[0] == r["decision"].split(":")[0], (wave, r["species"], r["decision"], new)
