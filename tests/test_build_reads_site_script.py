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


def test_screen_results_reach_the_species_report_and_the_counts(tmp_path):
    """Runs the worker screened (screen columns in the results table) carry them into the report; unscreened runs stay blank."""
    wave = make_repo(tmp_path, {}, ["MAT1-1"] * 6)
    f = wave / "results_all_final.tsv"
    rows = list(csv.DictReader(open(f, newline=""), delimiter="\t"))
    extra = {0: ("mixed", "Histoplasma", "0.34"), 1: ("genus_mismatch", "", ""), 2: ("ok", "", ""), 3: ("borderline", "Gigaspora", "0.0052")}
    for i, r in enumerate(rows):
        r["screen"], r["screen_second"], r["screen_explained"] = extra.get(i, ("", "", ""))
    tsv(f, list(rows[0]), [list(r.values()) for r in rows])
    d = build(tmp_path, tmp_path.parent)
    by = {r["run"]: r for r in d["rows"]}
    assert (by["SRR0"]["screen"], by["SRR0"]["screen_second"], by["SRR0"]["screen_explained"]) == ("mixed", "Histoplasma", "0.34")
    assert by["SRR4"]["screen"] == "" and by["SRR2"]["screen"] == "ok"
    assert d["screened"] == 4 and d["screen_flagged"] == 3          # ok is screened but not flagged
    idx = json.load(open(tmp_path / "out" / "index.json"))["species"][0]
    assert (idx["screened"], idx["screen_flagged"]) == (4, 3)
    assert "screen" in open(tmp_path / "out" / "sp_one.tsv").readline().split("\t")


def test_a_results_table_without_screen_columns_gives_zero_counts(tmp_path):
    make_repo(tmp_path, {}, ["MAT1-1"] * 6)
    d = build(tmp_path, tmp_path.parent)
    assert d["screened"] == 0 and d["screen_flagged"] == 0 and all(r["screen"] == "" for r in d["rows"])


# ---- the lineage of each species, for the grouped overview ----
import importlib.util


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(path.parent))
    spec.loader.exec_module(mod)
    return mod


def _genomes(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    tsv(path, ["genome", "species", "class_", "order", "family", "genus"], rows)


def test_a_species_is_placed_by_its_name_else_by_its_genus(tmp_path):
    brs = _load("build_reads_site", SCRIPT)
    _genomes(tmp_path / "results/asco/genomes.tsv", [
        ["g1", "Aspergillus flavus", "Eurotiomycetes", "Eurotiales", "Aspergillaceae", "Aspergillus"],
        ["g2", "Coccidioides immitis", "Eurotiomycetes", "Onygenales", "Onygenaceae", "Coccidioides"],
        ["g3", "Ascomycota sp.", "NA", "", "", "Ascomycota"],
    ])
    _genomes(tmp_path / "results/basidio/genomes.tsv", [["g4", "Cryptococcus neoformans", "Tremellomycetes", "Tremellales", "Cryptococcaceae", "Cryptococcus"]])
    table = brs.lineage_table(str(tmp_path), folders=[("asco", "Ascomycota"), ("basidio", "Basidiomycota"), ("missing", "Mucoromycota")])
    assert brs.lineage_fields(table, "Aspergillus flavus") == {"phylum": "Ascomycota", "class": "Eurotiomycetes", "order": "Eurotiales", "family": "Aspergillaceae"}
    # a name the tables lack lands with its genus, including one with a strain or collection note
    assert brs.lineage_fields(table, "Coccidioides posadasii")["order"] == "Onygenales"
    assert brs.lineage_fields(table, "Coccidioides immitis and C. posadasii")["family"] == "Onygenaceae"
    assert brs.lineage_fields(table, "Cryptococcus gattii (curated collection)")["phylum"] == "Basidiomycota"
    # an unknown genus is blank, and "NA" is a blank rank, not a class called NA
    assert brs.lineage_fields(table, "Nothing known") == {"phylum": "", "class": "", "order": "", "family": ""}
    assert brs.lineage_fields(table, "Ascomycota sp.") == {"phylum": "Ascomycota", "class": "", "order": "", "family": ""}


def test_a_missing_results_folder_gives_an_empty_table_and_the_build_still_runs(tmp_path):
    brs = _load("build_reads_site", SCRIPT)
    assert brs.lineage_table(str(tmp_path), folders=[("nowhere", "Ascomycota")]) == {}
    make_repo(tmp_path, {}, ["MAT1-1"] * 6)
    d = build(tmp_path, tmp_path.parent)
    assert (d["phylum"], d["class"], d["order"], d["family"]) == ("", "", "", "")


def test_the_species_report_and_the_index_carry_the_lineage(tmp_path):
    bsd = _load("build_site_data", REPO / "scripts" / "build_site_data.py")
    rel = bsd.RELEASES[bsd.DEFAULT_RELEASE]
    _genomes(tmp_path / "results" / rel["asco"] / "genomes.tsv", [["g1", "Species one", "Sordariomycetes", "Hypocreales", "Nectriaceae", "Species"]])
    make_repo(tmp_path, {}, ["MAT1-1"] * 6)
    d = build(tmp_path, tmp_path.parent)
    assert (d["phylum"], d["class"], d["order"], d["family"]) == ("Ascomycota", "Sordariomycetes", "Hypocreales", "Nectriaceae")
    idx = json.load(open(tmp_path / "out" / "index.json"))
    assert (idx["species"][0]["phylum"], idx["species"][0]["order"]) == ("Ascomycota", "Hypocreales")


def test_study_titles_reach_the_species_report_for_its_own_studies_only(tmp_path):
    """study_titles.tsv (from PopSeq) gives the hover text of the study links; a report carries only the titles of its own studies."""
    make_repo(tmp_path, {}, ["MAT1-1", "MAT1-2", "MAT1-1", "MAT1-2"])
    tsv(tmp_path / "datasets/reads/sra_survey/pilot/study_titles.tsv", ["study", "title"],
        [["PRJ1", "Genomic determinants of host specialization"], ["PRJ2", "Second study"], ["PRJ9", "A study of another species"]])
    d = build(tmp_path, tmp_path.parent)
    assert d["study_titles"] == {"PRJ1": "Genomic determinants of host specialization", "PRJ2": "Second study"}


def test_no_study_titles_file_gives_an_empty_map(tmp_path):
    make_repo(tmp_path, {}, ["MAT1-1", "MAT1-2"])
    assert build(tmp_path, tmp_path.parent)["study_titles"] == {}
