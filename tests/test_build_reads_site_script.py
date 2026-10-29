"""scripts/build_reads_site.py: per-species reports with clone groups, the plateau rule, and paths that do not depend on the working directory."""
import csv, json, subprocess, sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "build_reads_site.py"
RES = ["run", "group", "strain", "study", "layout", "state", "call", "truth", "verdict", "depth_1", "depth_2", "shared_depth", "breadth_1", "breadth_2", "flags", "truth_assembly", "idiomorphs"]


def tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(header); w.writerows(rows)


def make_repo(tmp, ani, calls):
    """A repo with one wave and one species of len(calls) runs; ani maps a run pair (by index) to its sketch identity (default 0.9)."""
    wave = tmp / "datasets/reads/sra_survey/pilot/waveX"; (wave / "clone_pairs").mkdir(parents=True)
    runs = [f"SRR{i}" for i in range(len(calls))]
    tsv(wave / "results_all_final.tsv", RES, [[r, "sp_one", f"s{i}", "PRJ1" if i % 2 else "PRJ2", "PAIRED", "done", c, "", "", "5", "0", "5", "0.9", "0", "-", "", "MAT1-1 MAT1-2"] for i, (r, c) in enumerate(zip(runs, calls))])
    tsv(wave / "gates_final.tsv", ["slug", "species", "decision", "reasons"], [["sp_one", "Species one", "PASS", "runs=%d" % len(runs)]])
    tsv(wave / "species.tsv", ["species", "slug", "strains_pooled", "reference_tier"], [["Species one", "sp_one", "40", "A"]])
    pairs = [[runs[i], runs[j], ani.get((i, j), 0.9)] for i in range(len(runs)) for j in range(i + 1, len(runs))]
    tsv(wave / "clone_pairs" / "sp_one.pairs.tsv", ["run_a", "run_b", "ani"], pairs)
    return wave


def build(tmp, cwd):
    out = tmp / "out"
    subprocess.run([sys.executable, str(SCRIPT), "--repo", str(tmp), "--wave", "datasets/reads/sra_survey/pilot/waveX", "--clones", "datasets/reads/sra_survey/pilot/waveX/clone_pairs", "--out", str(out)], check=True, cwd=cwd, capture_output=True)
    return json.load(open(out / "sp_one.json"))


def test_near_identical_runs_form_one_clone_and_plateau_gives_counts(tmp_path):
    calls = ["MAT1-1"] * 3 + ["MAT1-2"] * 9           # 12 runs; runs 0-2 are copies
    make_repo(tmp_path, {(0, 1): 0.9995, (0, 2): 0.9995, (1, 2): 0.9995}, calls)
    d = build(tmp_path, tmp_path.parent)              # a working directory other than the repo
    cl = d["clones"]
    assert cl["runs"] == 12 and cl["clones"] == 10 and len(cl["groups"]) == 1 and cl["groups"][0]["size"] == 3
    assert cl["plateau"] and cl["counts"] == {"MAT1-1": 1, "MAT1-2": 9}
    assert {r["clone"] for r in d["rows"] if r["run"] in ("SRR0", "SRR1", "SRR2")} == {"c1"}
    assert all(r["cid"] for r in d["rows"])


def test_unstable_clone_number_withholds_counts(tmp_path):
    # similar-but-not-identical pairs sit between 0.997 and 0.998: the clone number moves a lot between the two thresholds
    ani = {(i, i + 1): 0.9975 for i in range(0, 12, 2)}
    make_repo(tmp_path, ani, ["MAT1-1"] * 12)
    cl = build(tmp_path, tmp_path)["clones"]
    assert cl["clones"] == 12 and cl["clones_at_alt"] == 6 and cl["plateau"] is False and cl["counts"] is None


def test_few_runs_never_give_counts(tmp_path):
    make_repo(tmp_path, {(0, 1): 0.9995}, ["MAT1-1"] * 6)
    assert build(tmp_path, tmp_path)["clones"]["counts"] is None


def test_clone_mixing_single_and_both_calls_is_a_conflict(tmp_path):
    calls = ["MAT1-1", "both"] + ["MAT1-2"] * 10
    make_repo(tmp_path, {(0, 1): 0.9995}, calls)
    assert build(tmp_path, tmp_path)["clones"]["counts"].get("conflict") == 1
