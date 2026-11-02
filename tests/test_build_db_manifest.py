"""The browser page's database archive: complete, reproducible, and matching its manifest."""
import hashlib
import importlib.util
import io
import json
import sys
import tarfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("build_db_manifest", REPO / "web/spike/page/build_db_manifest.py")
bdm = importlib.util.module_from_spec(spec)
sys.modules["build_db_manifest"] = bdm
spec.loader.exec_module(bdm)


def _assets(tmp_path: Path, files: dict[str, bytes]) -> Path:
    assets = tmp_path / "_assets"
    for rel, data in files.items():
        p = assets / "db" / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)
    (assets / "bundle").mkdir(parents=True)
    (assets / "bundle" / "pyhmmer-0.11-cp312-none-any.whl").write_bytes(b"w")
    return assets


FILES = {"Ascomycota/order.yml": b"a: 1\n", "Ascomycota/Serinales/1_x/metadata.yaml": b"record_id: 1_x\n", "schema.json": b"{}"}


def test_the_archive_holds_every_database_file_under_its_relative_path(tmp_path):
    assets = _assets(tmp_path, FILES)
    m = bdm.write_assets(str(assets))
    with tarfile.open(assets / "db.tar.gz", "r:gz") as tar:
        got = {i.name: tar.extractfile(i).read() for i in tar.getmembers()}
    assert got == FILES
    assert m["db"] == sorted(FILES)


def test_the_manifest_matches_the_archive(tmp_path):
    assets = _assets(tmp_path, FILES)
    bdm.write_assets(str(assets))
    manifest = json.loads((assets / "db-manifest.json").read_text())
    blob = (assets / "db.tar.gz").read_bytes()
    assert manifest["dbArchive"] == "db.tar.gz"
    assert manifest["dbArchiveBytes"] == len(blob)
    assert manifest["dbArchiveSha256"] == hashlib.sha256(blob).hexdigest()
    assert manifest["pyhmmerWheel"] == "pyhmmer-0.11-cp312-none-any.whl"


def test_the_same_database_gives_the_same_bytes(tmp_path):
    a = bdm.write_assets(str(_assets(tmp_path / "one", FILES)))
    b = bdm.write_assets(str(_assets(tmp_path / "two", dict(reversed(list(FILES.items()))))))
    assert a["dbArchiveSha256"] == b["dbArchiveSha256"]


def test_no_timestamps_or_owners_leak_into_the_archive(tmp_path):
    assets = _assets(tmp_path, FILES)
    bdm.write_assets(str(assets))
    blob = (assets / "db.tar.gz").read_bytes()
    assert blob[4:8] == b"\x00\x00\x00\x00"  # the gzip header's mtime
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tar:
        assert {(i.mtime, i.uid, i.gid, i.uname) for i in tar.getmembers()} == {(0, 0, 0, "")}


def test_a_missing_pyhmmer_wheel_or_empty_database_stops_the_build(tmp_path):
    import pytest
    assets = _assets(tmp_path, FILES)
    (assets / "bundle" / "pyhmmer-0.11-cp312-none-any.whl").unlink()
    with pytest.raises(SystemExit, match="pyhmmer"):
        bdm.write_assets(str(assets))
    empty = tmp_path / "e" / "_assets"
    (empty / "db").mkdir(parents=True)
    (empty / "bundle").mkdir()
    (empty / "bundle" / "pyhmmer-1.whl").write_bytes(b"w")
    with pytest.raises(SystemExit, match="no database files"):
        bdm.write_assets(str(empty))


def test_the_real_database_archives_and_every_record_is_in_it(tmp_path):
    """The shipped db/ folder, copied as the site build does, round-trips through the archive."""
    import shutil
    assets = tmp_path / "_assets"
    shutil.copytree(REPO / "db", assets / "db")
    (assets / "bundle").mkdir()
    (assets / "bundle" / "pyhmmer-1.whl").write_bytes(b"w")
    m = bdm.write_assets(str(assets))
    with tarfile.open(assets / "db.tar.gz", "r:gz") as tar:
        assert sorted(tar.getnames()) == m["db"]
    assert len(m["db"]) == sum(1 for p in (REPO / "db").rglob("*") if p.is_file())
