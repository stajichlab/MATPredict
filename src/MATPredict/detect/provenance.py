"""The `run` block of `detection_report.yaml`: what was run, on what, with which code and database.

A report without it cannot say which genome it describes or whether two
reports used the same database (`curate-db release` does not stamp one yet),
so the per-genome HTML report (`matpredict report genome`) shows "not
recorded" for reports written before this block existed.
"""
from __future__ import annotations

import gzip
import hashlib
import time
from datetime import datetime, timezone
from pathlib import Path

from MATPredict import __version__


def genome_stats(path: Path) -> dict:
    """File name, SHA-256, contig count, total length and N50 of a FASTA (plain or .gz), in one pass."""
    digest = hashlib.sha256()
    lengths: list[int] = []
    current = None
    opener = gzip.open if str(path).endswith(".gz") else open
    with open(path, "rb") as raw:
        for block in iter(lambda: raw.read(1 << 20), b""):
            digest.update(block)
    with opener(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if current is not None:
                    lengths.append(current)
                current = 0
            elif current is not None:
                current += len(line.strip())
    if current is not None:
        lengths.append(current)
    total = sum(lengths)
    n50, acc = 0, 0
    for n in sorted(lengths, reverse=True):
        acc += n
        if acc * 2 >= total:
            n50 = n
            break
    return {"file": Path(path).name, "sha256": digest.hexdigest(), "contigs": len(lengths),
            "length_bp": total, "n50": n50}


def database_digest(db_root: Path) -> dict:
    """SHA-256 over the curated database files (path and content, sorted), and the record count.
    `db/candidates/` is left out: `detect` never reads it."""
    root = Path(db_root)
    digest = hashlib.sha256()
    records = 0
    for p in sorted(root.rglob("*")):
        rel = p.relative_to(root)
        if not p.is_file() or rel.parts[0] == "candidates":
            continue
        digest.update(str(rel).encode() + b"\0")
        digest.update(p.read_bytes())
        if p.name == "metadata.yaml":
            records += 1
    return {"root": str(root), "content_sha256": digest.hexdigest(), "records": records}


class RunClock:
    """Start time (UTC, ISO 8601) and elapsed wall seconds of one run."""

    def __init__(self) -> None:
        self.started = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
        self._t0 = time.monotonic()

    def wall_seconds(self) -> int:
        return int(round(time.monotonic() - self._t0))


def run_block(*, clock: RunClock, genome: Path | None, proteins: Path | None, db_root: Path,
              taxonomy_source: str, sample: str | None, organism: str | None, taxid: int | None,
              phylum: str | None, parameters: dict) -> dict:
    """The `run` mapping. A genome that cannot be read is reported, not raised: the report must still be written."""
    try:
        genome_doc = genome_stats(genome) if genome else None
    except OSError as exc:
        genome_doc = {"file": Path(genome).name, "error": f"{type(exc).__name__}: {exc}"}
    return {
        "sample": sample or (Path(genome).name.split(".")[0] if genome else None),
        "organism": organism,
        "taxid": taxid,
        "phylum": phylum,
        "matpredict_version": __version__,
        "database": database_digest(db_root),
        "taxonomy_source": taxonomy_source,
        "genome": genome_doc,
        "proteins": Path(proteins).name if proteins else None,
        "parameters": parameters,
        "started": clock.started,
        "wall_seconds": clock.wall_seconds(),
    }
