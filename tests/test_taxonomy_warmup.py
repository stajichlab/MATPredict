"""scripts/taxonomy_warmup.py: the panel scripts' taxonomy warm-up and snapshot record (#324, #342). Loaded by path (`scripts/` is not a package)."""
import importlib.util
import io
from pathlib import Path

import pytest

from MATPredict.db import local_taxonomy as lt

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "taxonomy_warmup.py"

# root - cellular - Eukaryota - Fungi - Ascomycota (phylum) - Pachysolen (genus) - P. tannophilus (code 26)
NODES = [(1, 1, "no rank", 1), (131567, 1, "cellular root", 1), (2759, 131567, "domain", 1), (4751, 2759, "kingdom", 1),
         (4890, 4751, "phylum", 1), (4917, 4890, "genus", 26), (4918, 4917, "species", 26)]
NAMES = {1: "root", 131567: "cellular organisms", 2759: "Eukaryota", 4751: "Fungi", 4890: "Ascomycota", 4917: "Pachysolen",
         4918: "Pachysolen tannophilus"}


def _dmp(rows):
    return "".join("\t|\t".join(str(x) for x in r) + "\t|\n" for r in rows)


def _load():
    spec = importlib.util.spec_from_file_location("taxonomy_warmup", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def offline_table(tmp_path, monkeypatch):
    d = tmp_path / "dump"
    d.mkdir()
    (d / "nodes.dmp").write_text(_dmp([(t, p, r, "", 0, 1, c, 1, 1, 1, 0, 0, "") for t, p, r, c in NODES]))
    (d / "names.dmp").write_text(_dmp([(t, n, "", "scientific name") for t, n in NAMES.items()]))
    (d / "merged.dmp").write_text(_dmp([(99999, 4918)]))
    out = tmp_path / "tax.tsv.zst"
    lt.build_table(d, out, snapshot="2026-10-01")
    monkeypatch.setenv("MATPREDICT_TAXONOMY", str(out))
    monkeypatch.setenv("MATPREDICT_OFFLINE", "1")
    monkeypatch.setenv("TAXONOMY_WARMUP_SLEEP", "0")
    return out


def test_taxids_in_the_snapshot_pass_including_a_merged_one(offline_table):
    err = io.StringIO()
    assert _load().check([4918, 99999], err=err) == []
    assert err.getvalue() == ""


def test_a_taxid_missing_from_the_snapshot_fails_loudly_offline(offline_table):
    err = io.StringIO()
    assert _load().check([4918, 123456, 5158], err=err) == [123456, 5158]
    assert "taxid 123456 failed" in err.getvalue() and "snapshot 2026-10-01" in err.getvalue()


def test_main_exit_code_and_the_list_of_failed_taxids(offline_table, monkeypatch, capsys):
    mod = _load()
    monkeypatch.setattr("sys.stdin", io.StringIO("4918\n5158\n\n"))
    assert mod.main([]) == 3
    assert "1 failed: 5158" in capsys.readouterr().err
    monkeypatch.setattr("sys.stdin", io.StringIO("4918\n"))
    assert mod.main([]) == 0


def test_describe_records_snapshot_path_and_hash(offline_table):
    text = _load().describe()
    assert "snapshot 2026-10-01" in text.splitlines()[0]
    assert "offline: True" in text and f"table: {offline_table}" in text and "snapshot: 2026-10-01" in text
    assert [l for l in text.splitlines() if l.startswith("sha256: ")][0].split()[1].isalnum()


def test_describe_says_so_without_a_table(monkeypatch):
    monkeypatch.delenv("MATPREDICT_TAXONOMY", raising=False)
    monkeypatch.setenv("MATPREDICT_OFFLINE", "1")
    assert "table: none" in _load().describe()
