"""Genomes the curator has suppressed are never run.

Curator's ruling 2026-09-26: the BFD suppress list
(Fungi_BFD/data/curation/suppress.txt, `ASMID,REASON`) applies here, plus a
MATPredict-only list in the curated database root. Nine PRJEB104476 records
are single Sanger amplicons deposited as "Complete Genome" (GCA_986280975.1 is
655 bp); one reached a panel as an empty FASTA and failed makeblastdb.
"""
import subprocess
import sys

from MATPredict.detect.suppress import filter_rows, load_suppress_list


def test_the_bfd_csv_format_is_read(tmp_path):
    f = tmp_path / "suppress.txt"
    f.write_text("ASMID,REASON\nGCA_986280975.1_DAH1005FM,broken upload\n"
                 "GCA_000091065.2_ASM9106v2,Too Small\n")
    assert load_suppress_list([f]) == {
        "GCA_986280975.1_DAH1005FM": "broken upload",
        "GCA_000091065.2_ASM9106v2": "Too Small",
    }


def test_one_id_per_line_with_comments_and_blanks(tmp_path):
    f = tmp_path / "local.txt"
    f.write_text("# MATPredict-only\n\nGCA_1.1_x\n  GCA_2.1_y  # trailing note\n")
    assert set(load_suppress_list([f])) == {"GCA_1.1_x", "GCA_2.1_y"}


def test_a_missing_file_is_skipped_not_fatal(tmp_path):
    assert load_suppress_list([tmp_path / "absent.txt"]) == {}


def test_lists_are_merged(tmp_path):
    a = tmp_path / "a.txt"; a.write_text("GCA_1.1_x\n")
    b = tmp_path / "b.txt"; b.write_text("GCA_2.1_y,reason\n")
    assert set(load_suppress_list([a, b])) == {"GCA_1.1_x", "GCA_2.1_y"}


def test_rows_are_matched_by_asmid_or_by_bare_accession():
    suppressed = {"GCA_986280975.1_DAH1005FM": "broken"}
    rows = ["GCA_986280975.1_DAH1005FM\t688394", "GCA_986280975.1\t688394",
            "GCA_000001.1_ok\t5"]
    kept, skipped = filter_rows(rows, suppressed)
    assert kept == ["GCA_000001.1_ok\t5"]
    assert len(skipped) == 2


def test_the_cli_filters_a_panel_list(tmp_path):
    sup = tmp_path / "suppress.txt"
    sup.write_text("ASMID,REASON\nGCA_9.1_bad,broken\n")
    panel = tmp_path / "panel.tsv"
    panel.write_text("GCA_9.1_bad\t1\nGCA_8.1_good\t2\n")
    proc = subprocess.run(
        [sys.executable, "-m", "MATPredict", "detect", "suppress-filter",
         "--list", str(panel), "--suppress", str(sup), "--no-default-lists"],
        capture_output=True, text=True, check=True,
    )
    assert proc.stdout == "GCA_8.1_good\t2\n"
    assert "skipped 1 suppressed" in proc.stderr


def test_the_allow_list_removes_genomes_from_the_suppress_set_by_full_or_bare_accession(tmp_path):
    """#165: rust genomes suppressed for annotation reasons are searched when they are on the MATPredict-only allow list."""
    bfd = tmp_path / "bfd.txt"
    bfd.write_text("ASMID,REASON\nGCA_902702905.1_Austropuccinia_psidii,EVM never completes\nGCA_943181505.1_UFV02.1,Too small (800kb)\nGCA_9.1_x,other\n")
    allow = tmp_path / "allow.txt"
    allow.write_text("# allowed\nGCA_902702905.1,Austropuccinia psidii: annotation only\n")
    got = load_suppress_list([bfd], allow=allow)
    assert set(got) == {"GCA_943181505.1_UFV02.1", "GCA_9.1_x"}                    # the allowed genome is gone, the legitimate skips stay
    kept, skipped = filter_rows(["GCA_902702905.1_Austropuccinia_psidii\t4", "GCA_943181505.1_UFV02.1\t7"], got)
    assert len(kept) == 1 and kept[0].startswith("GCA_902702905.1") and len(skipped) == 1


def test_a_missing_allow_list_changes_nothing(tmp_path):
    bfd = tmp_path / "bfd.txt"
    bfd.write_text("GCA_1.1_x,r\n")
    assert load_suppress_list([bfd], allow=tmp_path / "absent.txt") == load_suppress_list([bfd])


def test_the_committed_allow_list_names_only_genomes_that_are_on_the_bfd_reasons_and_not_the_too_small_one():
    from pathlib import Path
    ids = set(load_suppress_list([Path(__file__).resolve().parents[2] / "db" / "allow.txt"]))
    assert len(ids) == 9 and "GCA_943181505.1" not in ids and all(i.startswith("GCA_") for i in ids)
