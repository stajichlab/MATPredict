#!/usr/bin/env python3
"""Write the browser page's database archive and manifest into a built site's `_assets` folder.

The page used to request every database file on its own (767 requests per run), and one 503 from the
host in that burst ended the run. It now requests one archive, `db.tar.gz`, and checks it against the
manifest. The manifest still lists every file, so the page can count what it unpacked, and it names the
pyhmmer wheel.

Usage: build_db_manifest.py <assets-dir>   (needs <assets-dir>/db and <assets-dir>/bundle/pyhmmer-*.whl)

The archive is reproducible: sorted names, fixed times and owners, no gzip timestamp, so the same database
gives the same bytes and the same SHA-256.
"""
import gzip
import hashlib
import io
import json
import os
import sys
import tarfile

ARCHIVE = "db.tar.gz"


def database_files(db_dir: str) -> list[str]:
    """Relative paths of every file under `db_dir`, sorted, with `/` separators."""
    out = []
    for root, _, names in os.walk(db_dir):
        for n in names:
            out.append(os.path.relpath(os.path.join(root, n), db_dir).replace(os.sep, "/"))
    return sorted(out)


def build_archive(db_dir: str, files: list[str]) -> bytes:
    """A gzipped tar of `files`, relative to `db_dir`, the same bytes every time."""
    raw = io.BytesIO()
    with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, compresslevel=9) as gz:
        with tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as tar:
            for rel in files:
                path = os.path.join(db_dir, *rel.split("/"))
                info = tarfile.TarInfo(rel)
                info.size = os.path.getsize(path)
                info.mtime = 0
                info.mode = 0o644
                info.uid = info.gid = 0
                info.uname = info.gname = ""
                with open(path, "rb") as fh:
                    tar.addfile(info, fh)
    return raw.getvalue()


def write_assets(assets: str) -> dict:
    """Write `db.tar.gz` and `db-manifest.json` under `assets`; return the manifest."""
    db_dir = os.path.join(assets, "db")
    bundle = os.path.join(assets, "bundle")
    files = database_files(db_dir)
    wheels = sorted(n for n in os.listdir(bundle) if n.startswith("pyhmmer-") and n.endswith(".whl"))
    if not wheels:
        raise SystemExit("no pyhmmer wheel in the bundle")
    if not files:
        raise SystemExit(f"no database files under {db_dir}")
    blob = build_archive(db_dir, files)
    with open(os.path.join(assets, ARCHIVE), "wb") as fh:
        fh.write(blob)
    manifest = {
        "db": files,
        "pyhmmerWheel": wheels[0],
        "dbArchive": ARCHIVE,
        "dbArchiveBytes": len(blob),
        "dbArchiveSha256": hashlib.sha256(blob).hexdigest(),
    }
    with open(os.path.join(assets, "db-manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
    return manifest


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    m = write_assets(argv[1])
    print(f"manifest: {len(m['db'])} database files, {m['pyhmmerWheel']}; "
          f"{m['dbArchive']} {m['dbArchiveBytes'] / 1024:.0f} KB sha256 {m['dbArchiveSha256'][:12]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
