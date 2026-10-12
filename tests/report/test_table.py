"""`matpredict report table`: one row per genome and one row per locus, from a directory of detect runs."""
from __future__ import annotations

import csv
import shutil
from pathlib import Path

from MATPredict.__main__ import main

FIXTURES = Path(__file__).parent / "fixtures"


def _runs(tmp_path: Path, **name_to_fixture: str) -> Path:
    """A runs directory: one run folder per name, holding that fixture as detection_report.yaml."""
    runs = tmp_path / "runs"
    for name, fixture in name_to_fixture.items():
        (runs / name).mkdir(parents=True)
        shutil.copy(FIXTURES / f"{fixture}.yaml", runs / name / "detection_report.yaml")
    return runs


def _tsv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_genomes_table_has_one_row_per_run(tmp_path):
    runs = _runs(tmp_path, g_two="two_idiomorphs", g_none="none_called_gap", g_skip="not_searched")
    out = tmp_path / "out"
    assert main(["report", "table", "--runs", str(runs), "--out-dir", str(out)]) == 0
    rows = {r["sample"]: r for r in _tsv(out / "genomes.tsv")}
    assert set(rows) == {"g_two", "g_none", "g_skip"}
    assert rows["g_two"]["status"] == "report"
    assert rows["g_two"]["n_loci"] == "2"
    assert rows["g_two"]["idiomorphs"] == "Minus;Plus"
    assert rows["g_none"]["n_loci"] == "0"
    assert rows["g_none"]["n_withheld"] == "2"
    assert rows["g_skip"]["routing_mode"] == "not_searched"
    assert "Chytridiomycota" in rows["g_skip"]["not_searched_reason"]


def test_loci_table_has_one_row_per_called_locus(tmp_path):
    runs = _runs(tmp_path, g_two="two_idiomorphs", g_none="none_called_gap")
    out = tmp_path / "out"
    assert main(["report", "table", "--runs", str(runs), "--out-dir", str(out)]) == 0
    loci = _tsv(out / "loci.tsv")
    assert [r["sample"] for r in loci] == ["g_two", "g_two"]
    plus = next(r for r in loci if r["idiomorph"] == "Plus")
    assert (plus["family"], plus["contig"], plus["start"], plus["end"]) == ("Mucoromycota:MAT", "scaffold_12", "50100", "56300")
    assert (plus["confidence"], plus["locus_class"]) == ("high", "mat_locus")
    minus = next(r for r in loci if r["idiomorph"] == "Minus")
    assert (minus["contig"], minus["start"], minus["end"], minus["confidence"]) == ("scaffold_88", "1200", "1980", "low")


def test_runs_without_a_report_still_get_a_row_with_the_reason(tmp_path):
    runs = _runs(tmp_path, g_ok="none_called_gap")
    (runs / "g_fail").mkdir()
    (runs / "g_fail" / "failed_reason.txt").write_text("rc=124: timed out\n")
    (runs / "g_skip").mkdir()
    (runs / "g_skip" / "skipped_reason.txt").write_text("taxid not in the table\n")
    (runs / "g_empty").mkdir()
    out = tmp_path / "out"
    assert main(["report", "table", "--runs", str(runs), "--out-dir", str(out)]) == 0
    rows = {r["sample"]: r for r in _tsv(out / "genomes.tsv")}
    assert (rows["g_fail"]["status"], rows["g_fail"]["reason"]) == ("failed", "rc=124: timed out")
    assert (rows["g_skip"]["status"], rows["g_skip"]["reason"]) == ("skipped", "taxid not in the table")
    assert rows["g_empty"]["status"] == "no_report"
    assert rows["g_ok"]["status"] == "report"
    assert rows["g_ok"]["reason"] == ""


def test_genome_row_carries_assembly_run_and_mating_type_signals(tmp_path):
    runs = _runs(tmp_path, g_two="two_idiomorphs", g_none="none_called_gap")
    out = tmp_path / "out"
    assert main(["report", "table", "--runs", str(runs), "--out-dir", str(out)]) == 0
    rows = {r["sample"]: r for r in _tsv(out / "genomes.tsv")}
    two = rows["g_two"]
    assert (two["organism"], two["taxid"], two["phylum"]) == ("Rhizopus stolonifer", "4846", "Mucoromycota")
    assert (two["contigs"], two["length_bp"], two["n50"], two["wall_seconds"]) == ("812", "48213377", "1203321", "96")
    assert two["genetic_code"] == "1"
    assert two["families_attempted"] == "Mucoromycota:MAT"
    # the report's own two_idiomorphs record (arrangement), and the one-locus homothallic class
    assert two["two_idiomorphs"] == "unlinked"
    assert two["homothallic_candidate"] == "0"
    assert rows["g_none"]["two_idiomorphs"] == ""


def test_samples_file_adds_metadata_and_rows_for_genomes_never_run(tmp_path):
    runs = _runs(tmp_path, g_two="two_idiomorphs")
    samples = tmp_path / "samples.csv"
    samples.write_text("ASMID,ORDER,SPECIES\ng_two,Mucorales,Rhizopus stolonifer\ng_absent,Pezizales,Morchella sp.\n")
    out = tmp_path / "out"
    assert main(["report", "table", "--runs", str(runs), "--samples", str(samples), "--out-dir", str(out)]) == 0
    rows = {r["sample"]: r for r in _tsv(out / "genomes.tsv")}
    assert rows["g_two"]["ORDER"] == "Mucorales"
    assert (rows["g_absent"]["status"], rows["g_absent"]["ORDER"]) == ("not_run", "Pezizales")
    assert _tsv(out / "loci.tsv")[0]["sample"] == "g_two"


def _table_and_summary(tmp_path, by="PHYLUM"):
    runs = _runs(tmp_path, g_two="two_idiomorphs", g_none="none_called_gap", g_skip="not_searched", g_ok="basidio_hd_pr")
    (runs / "g_fail").mkdir()
    (runs / "g_fail" / "failed_reason.txt").write_text("rc=124\n")
    samples = tmp_path / "samples.csv"
    samples.write_text("ASMID,PHYLUM\ng_two,Mucoromycota\ng_none,Ascomycota\ng_skip,Chytridiomycota\ng_ok,Basidiomycota\ng_fail,Ascomycota\ng_absent,Ascomycota\n")
    out = tmp_path / "out"
    assert main(["report", "table", "--runs", str(runs), "--samples", str(samples), "--out-dir", str(out)]) == 0
    assert main(["report", "summary", "--table-dir", str(out), "--by", by]) == 0
    return out


def test_summary_counts_recovery_by_group(tmp_path):
    out = _table_and_summary(tmp_path)
    rows = {r["PHYLUM"]: r for r in _tsv(out / "summary_by_PHYLUM.tsv")}
    asco = rows["Ascomycota"]
    # g_none searched, no call; g_fail failed; g_absent never run
    assert (asco["genomes"], asco["searched"], asco["called"], asco["failed"], asco["not_run"]) == ("3", "1", "0", "1", "1")
    muc = rows["Mucoromycota"]
    assert (muc["searched"], muc["called"], muc["recovery"]) == ("1", "1", "1.00")
    # a phylum with no curated family is counted as not searched, not as a miss
    chy = rows["Chytridiomycota"]
    assert (chy["searched"], chy["not_searched"], chy["recovery"]) == ("0", "1", "")
    assert "Mucoromycota" in (out / "summary_by_PHYLUM.md").read_text()


def test_summary_counts_idiomorphs_and_both_idiomorph_genomes(tmp_path):
    out = _table_and_summary(tmp_path)
    rows = {r["PHYLUM"]: r for r in _tsv(out / "summary_by_PHYLUM.tsv")}
    muc = rows["Mucoromycota"]
    assert set(muc["idiomorph_loci"].split(";")) == {"Plus:1", "Minus:1"}
    assert muc["both_idiomorphs"] == "1"
    assert muc["called_high_confidence"] == "1"
    allrow = rows["ALL"]
    assert (allrow["genomes"], allrow["not_run"], allrow["failed"], allrow["not_searched"]) == ("6", "1", "1", "1")


def test_summary_by_a_missing_column_says_which_columns_exist(tmp_path, caplog):
    runs = _runs(tmp_path, g_two="two_idiomorphs")
    out = tmp_path / "out"
    main(["report", "table", "--runs", str(runs), "--out-dir", str(out)])
    assert main(["report", "summary", "--table-dir", str(out), "--by", "ORDER"]) != 0
    assert "no column 'ORDER'" in caplog.text
