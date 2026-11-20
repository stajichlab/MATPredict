"""Fixture test of scripts/model_disagreement_summary.py (#290, #349)."""
from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("mds", Path(__file__).resolve().parents[1] / "scripts" / "model_disagreement_summary.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

HEAD = ["class", "role", "trigger_noncanonical", "primary_identity", "alternate_identity", "primary_aa", "alternate_aa", "ref_aa"]


def write(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(HEAD)
        w.writerows(rows)


def test_yardsticks_counts_and_medians():
    rows = [{"alternate_identity": "90", "primary_identity": "80", "alternate_aa": "100", "primary_aa": "80", "ref_aa": "100"},
            {"alternate_identity": "70", "primary_identity": "70", "alternate_aa": "90", "primary_aa": "100", "ref_aa": "100"},
            {"alternate_identity": "60", "primary_identity": "75", "alternate_aa": "50", "primary_aa": "100", "ref_aa": "100"},
            {"alternate_identity": "", "primary_identity": "75", "alternate_aa": "50", "primary_aa": "50", "ref_aa": ""}]
    y = mod.yardsticks(rows)
    assert y["genes"] == 4 and y["with_identity"] == 3 and y["alt_higher_identity"] == 1 and y["identity_tie"] == 1
    assert y["median_identity_diff"] == 0 and y["with_ref_length"] == 3 and y["alt_closer_length"] == 1 and y["length_tie"] == 0 and y["median_aa_diff"] == -10.0


def test_run_subsets_and_skips_not_checked(tmp_path):
    write(tmp_path / "x.census.tsv", [["primary_issue_alt_clean", "core_MAT", "True", "70", "90", "80", "100", "100"],
                                      ["both_clean", "flanking_conserved", "False", "80", "80", "100", "100", "100"],
                                      ["not_checked", "core_MAT", "False", "", "", "", "", ""]])
    rows = mod.run([tmp_path / "x.census.tsv"], tmp_path / "o.tsv")
    by = {r["subset"]: r for r in rows}
    assert by["all"]["genes"] == 2 and by["core_MAT"]["genes"] == 1 and by["trigger_noncanonical"]["genes"] == 1
    assert by["primary_issue_alt_clean"]["alt_higher_identity"] == 1 and by["primary_issue_alt_clean|core_MAT"]["genes"] == 1
    assert (tmp_path / "o.tsv").read_text().splitlines()[0].startswith("run\tsubset\tgenes")
