"""reads_availability.py: PopSeq's run table joined to typed runs, gates and local data; snapshot and CHANGES.md."""
import csv
import gzip
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "datasets/reads/availability/reads_availability.py"
sys.path.insert(0, str(ROOT / "datasets/reads/sra_survey"))
from survey_ena import classify_study  # noqa: E402

RUN_COLS = ["run_accession", "study_accession", "sample_accession", "scientific_name", "tax_id", "library_layout", "read_count",
            "base_count", "strain", "isolate", "sample_title", "study_title"]


def tsv(p, cols, rows, delim="\t"):
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter=delim, lineterminator="\n")
        w.writerow(cols)
        w.writerows(rows)


def test_classify_study_keeps_the_survey_rules():
    assert classify_study("Genome sequencing of knockout strains", {"s1": "x"}) == "mutant"
    assert classify_study("Population genomics of clinical isolates", {"s1": "progeny 1"}) == ""     # a protected title is not 'cross' on one strain word
    assert classify_study("A study", {"a": "F1 progeny", "b": "segregant 2"}) == "cross"
    assert classify_study("RNA-seq time course", {"a": "x"}) == "not_wgs"
    assert classify_study("", {}) == ""


def fixture(tmp):
    big, small = "3000000000", "100000000"
    runs = [
        # Fungus alpha (taxid 11): natural study with 4 BioSamples (B1 has two runs), and a mutant study
        ["R1", "S1", "B1", "Fungus alpha", "11", "PAIRED", "10", big, "st1", "", "", "Natural isolates of F. alpha"],
        ["R1b", "S1", "B1", "Fungus alpha", "11", "PAIRED", "10", "20000000000", "st1", "", "", "Natural isolates of F. alpha"],
        ["R2", "S1", "B2", "Fungus alpha", "11", "PAIRED", "10", big, "st2", "", "", "Natural isolates of F. alpha"],
        ["R3", "S1", "B3", "Fungus alpha", "11", "SINGLE", "10", big, "st3", "", "", "Natural isolates of F. alpha"],
        ["R4", "S1", "B4", "Fungus alpha var. beta", "12", "PAIRED", "10", big, "st4", "", "", "Natural isolates of F. alpha"],
        ["R5", "S2", "B5", "Fungus alpha", "11", "PAIRED", "10", big, "m1", "", "", "Knockout collection"],
        ["R6", "S1", "B6", "Fungus alpha", "11", "PAIRED", "10", small, "st6", "", "", "Natural isolates of F. alpha"],   # below 300 Mb
        # Fungus gamma (taxid 21): 5 isolates, never typed
        *[[f"G{i}", "S3", f"GB{i}", "Fungus gamma", "21", "PAIRED", "10", big, f"g{i}", "", "", "Survey"] for i in range(5)],
    ]
    with gzip.open(tmp / "fungi_runs_ena_2030-01-01.tsv.gz", "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(RUN_COLS); w.writerows(runs)
    tsv(tmp / "taxmap.csv", ["taxid", "name", "rank", "species_taxid", "species_name", "resolution"],
        [["11", "Fungus alpha", "species", "11", "Fungus alpha", "self"], ["12", "Fungus alpha var. beta", "varietas", "11", "Fungus alpha", "rolled_up"],
         ["21", "Fungus gamma", "species", "21", "Fungus gamma", "self"]], delim=",")
    tsv(tmp / "typed.tsv", ["run", "group", "state", "call"], [["R1b", "fungus_alpha", "done", "MAT1-1"], ["R2", "fungus_alpha", "done", "MAT1-2"]])
    tsv(tmp / "gates.tsv", ["slug", "species", "decision", "reasons"], [["fungus_alpha", "Fungus alpha", "PASS", "ok"]])
    tsv(tmp / "cocci.tsv", ["strain", "species", "reads_call"], [])
    inv = tmp / "inv" / "panels"
    for sub, rows in (("Alpha_pop", [["st3", "", "", "local_only", "2"], ["new9", "", "", "local_only", "2"], ["gone", "", "", "listed_absent", "0"]]),
                      ("Alpha_lab", [["lab1", "", "", "local_only", "2"]])):
        (inv / sub).mkdir(parents=True)
        tsv(inv / sub / "metadata.csv", ["sample", "biosample", "run_accession", "provenance", "n_files"], rows, delim=",")
    tsv(tmp / "map.tsv", ["subproject", "species", "population", "note"],
        [["Alpha_pop", "Fungus alpha", "yes", ""], ["Alpha_lab", "Fungus alpha", "no", "lab lineage"]])


def run(tmp, *extra):
    cmd = [sys.executable, str(SCRIPT), "--runs", str(tmp / "fungi_runs_ena_2030-01-01.tsv.gz"), "--taxmap", str(tmp / "taxmap.csv"),
           "--typed", str(tmp / "typed.tsv"), "--gates", str(tmp / "gates.tsv"), "--typed-local", str(tmp / "cocci.tsv"),
           "--local-inventory", str(tmp / "inv"), "--local-map", str(tmp / "map.tsv"), "--out", str(tmp / "snap"), *extra]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    out = tmp / "snap" / "2030-01-01"
    return {r["species"]: r for r in csv.DictReader(open(out / "species_availability.tsv"), delimiter="\t")}, out


def test_species_rows_untyped_runs_and_local_counts(tmp_path):
    fixture(tmp_path)
    rows, out = run(tmp_path)
    a, g = rows["Fungus alpha"], rows["Fungus gamma"]
    # B1-B4 eligible (B4 rolled up from the variety; B5 mutant study; B6 too small); B1 and B2 typed
    assert (a["public_isolates_eligible"], a["typed_public"], a["untyped_public_eligible"]) == ("4", "2", "2")
    assert a["excluded_isolates"] == "mutant:1" and a["action"] == "expand" and a["gate"] == "PASS"
    # local: st3 names a public isolate, new9 does not, 'gone' has no files, the lab-lineage set is not counted
    assert (a["local_with_files"], a["local_likely_public"], a["local_not_public"]) == ("2", "1", "1")
    assert a["calls"] == "MAT1-1:1;MAT1-2:1" and a["minor_idiomorph_frac"] == "0.500"
    assert g["action"] == "new_candidate" and g["untyped_public_eligible"] == "5"      # no --db: tier unknown
    with gzip.open(out / "untyped_runs.tsv.gz", "rt") as fh:
        ut = list(csv.DictReader(fh, delimiter="\t"))
    assert sorted(r["run_accession"] for r in ut) == ["R3", "R4"]                    # only the expand species' untyped isolates
    local = list(csv.DictReader(open(out / "local_samples.tsv"), delimiter="\t"))
    assert {r["sample"]: r["counted"] for r in local}["lab1"] == "no"


def test_a_rerun_is_identical_and_changes_md_reports_growth(tmp_path):
    fixture(tmp_path)
    _, out = run(tmp_path)
    first = {p.name: p.read_bytes() for p in out.iterdir()}
    prev = tmp_path / "prev"
    out.rename(prev)
    # shrink the previous snapshot's count for gamma so the new run shows growth
    rows = list(csv.DictReader(open(prev / "species_availability.tsv"), delimiter="\t"))
    for r in rows:
        if r["species"] == "Fungus gamma":
            r["public_isolates_eligible"] = "3"
    with open(prev / "species_availability.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    _, out = run(tmp_path, "--previous", str(prev))
    assert {p.name: p.read_bytes() for p in out.iterdir() if p.name != "CHANGES.md"} == first
    assert "*Fungus gamma* (new_candidate): +2" in (out / "CHANGES.md").read_text()
