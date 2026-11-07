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
            new, _ = wd.decide_species(r["species"], counted, done, typed, both, truth, agree, msd)
            assert new.split(":")[0] == r["decision"].split(":")[0], (wave, r["species"], r["decision"], new)


def test_submit_full_takes_the_species_that_are_not_held_and_can_be_limited_to_one_decision(tmp_path, monkeypatch):
    import argparse
    wave = tmp_path / "wave"
    for sl, decision in (("p", "PASS"), ("u", "UNVALIDATED:typed_0.70"), ("h", "HOLD:typed_0.40"), ("v", "PROVISIONAL")):
        (wave / "plan" / sl).mkdir(parents=True)
        with open(wave / "plan" / sl / "runs_full.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["run_accession", "group"]); w.writerows([[f"{sl}{i}", sl] for i in range(5)])
    with open(wave / "gates.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["slug", "species", "decision", "reasons"])
        w.writerows([[s, s, d, ""] for s, d in (("p", "PASS"), ("u", "UNVALIDATED:typed_0.70"), ("h", "HOLD:typed_0.40"), ("v", "PROVISIONAL"))])
    seen = []
    monkeypatch.setattr(wd, "submit_stage", lambda a, stage, rows: seen.append((stage, [r["group"] for r in rows])))
    base = dict(wave=str(wave), cap=3, only=None, stage_name="list", gates="gates.tsv")
    wd.submit_full(argparse.Namespace(**base))
    assert seen[-1][0] == "full" and sorted(set(seen[-1][1])) == ["p", "u", "v"] and len(seen[-1][1]) == 9      # held species left out; at most --cap runs each
    wd.submit_full(argparse.Namespace(**{**base, "only": "UNVALIDATED", "stage_name": "full_unvalidated"}))
    assert seen[-1] == ("full_unvalidated", ["u", "u", "u"])
    import pytest
    with pytest.raises(SystemExit, match="no species to run with the decision NOT_A_DECISION"):
        wd.submit_full(argparse.Namespace(**{**base, "only": "NOT_A_DECISION"}))


def _res(call, truth="", state="done", sd="15"):
    return {"state": state, "call": call, "truth": truth, "shared_depth": sd}


def test_zygosity_species_count_both_as_a_call_and_as_agreeing_with_a_collapsed_assembly():
    # Candida albicans test stage of wave2: 7 single calls, 8 both, 5 low/none; with a truth: 2 both (truth a), 2 agree, 3 disagree, 2 low_depth
    rs = [_res("both", "a"), _res("both", "a"), _res("alpha", "alpha"), _res("alpha", "alpha"), _res("a", "alpha"), _res("a", "alpha"), _res("alpha", "a")] \
        + [_res("both")] * 6 + [_res("alpha"), _res("a")] + [_res("low_depth", "a")] * 3 + [_res("none", "alpha"), _res("low_depth", "alpha")]
    assert wd.numbers(rs)[2:6] == (7, 8, 5, 2)                                     # the ordinary count: typed 7, both 8, 5 truth runs, 2 agree
    n = wd.numbers(rs, zygosity=True)
    assert n[2:6] == (15, 0, 7, 4)                                                 # typed = single + both; both agrees with a single truth
    assert wd.decide_species("Candida albicans", *n)[0] == "UNVALIDATED:zygosity+agree_4of7"
    assert wd.decide_species("Candida albicans", 20, 20, 19, 0, 10, 10, 15.0)[0] == "PASS_ZYGOSITY"
    assert wd.decide_species("Candida albicans", 20, 20, 9, 0, 3, 3, 15.0)[0].startswith("HOLD:typed_0.45+provisional_typed_0.45+zygosity")
    assert wd.decide_species("Candida albicans", 20, 20, 2, 18, 0, 0, 15.0)[0] != "PROVISIONAL_HOMOTHALLIC"      # never the homothallic shortcut


def test_released_species_run_flagged_unless_the_downloads_failed():
    assert wd.decide_species("Saccharomyces paradoxus", 20, 20, 0, 0, 0, 0, 1.5)[0] == "UNVALIDATED:released_cassettes+typed_0.00+provisional_typed_0.00+shared_depth_1.5"
    assert wd.decide_species("Saccharomyces paradoxus", 20, 10, 0, 0, 0, 0, 1.5)[0].startswith("HOLD:many_failed")
    assert wd.decide_species("Species x", 20, 20, 0, 0, 0, 0, 1.5)[0].startswith("HOLD")                              # not on the list: unchanged


def test_curate_applies_the_rulings_and_leaves_other_species_alone(tmp_path):
    import argparse
    d = _gates(tmp_path, ROWS + [
        ["saccharomyces_paradoxus", "Saccharomyces paradoxus", "HOLD:typed_0.00+shared_depth_1.5+provisional_typed_0.00", "runs=20 done=20 typed=0 both=0 none/low=20 truth_typed=0 agree=0 median_shared_depth=1.5"],
        ["candida_albicans", "Candida albicans", "HOLD:typed_0.35+both_0.40+agree_2of5", "runs=4 done=4 typed=1 both=3 none/low=0 truth_typed=1 agree=0 median_shared_depth=15.9"]])
    (d / "test" / "summary").mkdir(parents=True)
    with open(d / "test" / "summary" / "results_all.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["run", "group", "state", "call", "truth", "shared_depth"])
        w.writerows([["r1", "candida_albicans", "done", "both", "a", "15"], ["r2", "candida_albicans", "done", "both", "", "15"], ["r3", "candida_albicans", "done", "both", "", "15"], ["r4", "candida_albicans", "done", "alpha", "", "15"]])
    before = (d / "gates_final.tsv").read_text().splitlines()
    wd.curate(argparse.Namespace(wave=str(d), in_file="gates_final.tsv", out_file=None, stage="test"))
    after = (d / "gates_final.tsv").read_text().splitlines()
    assert before[:len(ROWS) + 1] == after[:len(ROWS) + 1]
    assert after[-2].split("\t")[2].startswith("UNVALIDATED:released_cassettes")
    assert after[-1].split("\t")[2] == "PROVISIONAL_ZYGOSITY" and after[-1].endswith("mode=zygosity")
    _run_regate(d)                                                                 # regate keeps both rulings
    assert (d / "gates_final.tsv").read_text().splitlines()[-2:] == after[-2:]


def test_incomplete_lists_killed_missing_and_failed_runs_but_not_typing_errors(tmp_path):
    import argparse, json as js
    w = tmp_path / "wave"; (w / "st").mkdir(parents=True); (w / "out").mkdir()
    states = {"r1": "done", "r2": "started", "r3": "failed_download", "r4": "failed_type", "r5": "unavailable_at_ncbi"}      # r6: no run.json
    with open(w / "st" / "runs_all.tsv", "w", newline="") as fh:
        cw = csv.writer(fh, delimiter="\t", lineterminator="\n"); cw.writerow(["run_accession", "group"]); cw.writerows([[r, "g"] for r in list(states) + ["r6"]])
    for r, s in states.items():
        (w / "out" / f"{r}.run.json").write_text(js.dumps({"state": s}))
    wd.incomplete(argparse.Namespace(wave=str(w), stage_name="st", out_dir=None))
    got = [r["run_accession"] for r in csv.DictReader(open(w / "st" / "runs_incomplete.tsv"), delimiter="\t")]
    assert got == ["r2", "r3", "r6"]


def test_retry_pass_and_summary_never_go_to_the_run_partitions(tmp_path, monkeypatch):
    import argparse
    calls = []
    monkeypatch.setattr(wd, "sbatch", lambda args, **kw: calls.append(args) or str(len(calls)))
    a = argparse.Namespace(wave=str(tmp_path / "w"), out_dir=None, parts="preempt:preempt,short:stajichlab:short", retry_part="stajichlab:stajichlab", cpus=1, mem="8gb", conc=60, env=[], dry_run=False)
    wd.submit_stage(a, "st", [{"run_accession": f"r{i}", "group": "g"} for i in range(4)])
    part = lambda c: (c[c.index("-p") + 1], c[c.index("--qos") + 1] if "--qos" in c else "")
    names = [c[c.index("-J") + 1] for c in calls]
    assert [part(c) for c, n in zip(calls, names) if n.endswith("-run")] == [("preempt", ""), ("short", "short")]
    assert all(part(c)[0] == "stajichlab" for c, n in zip(calls, names) if n.endswith("-retry") or n.endswith("-summary"))
