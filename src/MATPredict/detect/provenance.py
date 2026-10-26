"""The `run` block of `detection_report.yaml`: what was run, on what, with which code and database.

A report without it cannot say which genome it describes or whether two
reports used the same database (`curate-db release` does not stamp one yet),
so the per-genome HTML report (`matpredict report genome`) shows "not
recorded" for reports written before this block existed.
"""
from __future__ import annotations

import functools
import gzip
import hashlib
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

from MATPredict import __version__, logger

try:  # the project requires Python >= 3.11; the guard keeps a bare import working on older interpreters
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    tomllib = None


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


def _git(top: Path, *args: str) -> str | None:
    try:
        done = subprocess.run(["git", "-C", str(top), *args], capture_output=True, text=True, errors="replace",
                              timeout=15, check=False)
    except (OSError, subprocess.SubprocessError, ValueError):
        return None
    return done.stdout.strip() if done.returncode == 0 else None


#: What counts as "the code" for `dirty`: the package, the curated database, the scripts the campaigns run, and the
#: project file. Untracked files count (a module dropped into src/ changes behaviour); results written elsewhere in
#: the tree, and docs, do not, so a tree does not turn dirty because a results file was updated mid-campaign.
CODE_PATHS = ("src", "db", "scripts", "pyproject.toml")


def source_version_for(package_dir: Path) -> str | None:
    """The `version` in the `pyproject.toml` of the source checkout that `package_dir` (the `MATPredict` package
    directory) sits in, or None: a wheel install (no `src/` layout), no `pyproject.toml`, another project's, or an
    unreadable one. Needs no git."""
    package_dir = Path(package_dir).resolve()
    if tomllib is None or package_dir.parent.name != "src":
        return None
    try:
        project = tomllib.loads((package_dir.parents[1] / "pyproject.toml").read_text(encoding="utf-8"))["project"]
        if str(project.get("name", "")).lower() != "matpredict":
            return None
        return str(project["version"])
    except (OSError, KeyError, ValueError):
        return None


@functools.lru_cache(maxsize=None)
def running_version() -> str:
    """The version a report records as `matpredict_version`: that of the source checkout the code runs from when there
    is one, else the installed package metadata. The metadata alone is stale whenever an editable install was not
    refreshed (issue #139: reports said 0.6.1 for code at 0.6.3), and a shared environment must not be reinstalled to fix
    it, so the checkout's own `pyproject.toml` is the better witness."""
    return source_version_for(Path(__file__).resolve().parents[1]) or __version__


def code_stamp_for(package_dir: Path) -> dict:
    """Identity of the source checkout `package_dir` (the installed `MATPredict` package directory) runs from.

    The installed package metadata is stale whenever the editable install was not refreshed (issue #139: the pilot's
    reports said 0.6.1 for code at 0.6.3). This reads the checkout itself: `git_describe` (`git describe --tags --always`: the nearest tag), the full `commit`, whether
    the code is `dirty` (any change or untracked file under `CODE_PATHS`), and the `source_version` in its
    `pyproject.toml`. Only keys that could be read are returned;
    `{}` when the package is not in a git checkout of this project (a wheel install, no git, an unrelated
    enclosing repository), so a report from such an install is unchanged.
    """
    package_dir = Path(package_dir).resolve()
    top_text = _git(package_dir, "rev-parse", "--show-toplevel")
    if not top_text:
        return {}
    top = Path(top_text)
    if (top / "src" / "MATPredict").resolve() != package_dir:
        return {}  # a repository around somebody's venv, not this project's checkout
    out: dict = {}
    describe = _git(top, "describe", "--tags", "--always")
    if describe:
        out["git_describe"] = describe
    changed = _git(top, "status", "--porcelain", "--untracked-files=normal", "--", *CODE_PATHS)
    if changed is not None:
        out["dirty"] = bool(changed)
    commit = _git(top, "rev-parse", "HEAD")
    if commit:
        out["commit"] = commit
    version = source_version_for(package_dir)
    if version:
        out["source_version"] = version
    return out


@functools.lru_cache(maxsize=None)
def code_stamp() -> dict:
    """`code_stamp_for` the running package, computed once per process. `detect` calls it at the start of a run, so
    the stamp describes the code that was checked out when the run began. A source checkout whose stamp cannot be read
    (git missing, a `safe.directory` refusal on a shared directory, a timeout) is logged, not silent: the report would
    otherwise lack `code` without anyone noticing."""
    package_dir = Path(__file__).resolve().parents[1]
    stamp = code_stamp_for(package_dir)
    if not stamp and (package_dir.parents[1] / ".git").exists():
        logger.warning("could not read the git identity of the source checkout %s; the report will carry no `code` stamp",
                       package_dir.parents[1])
    return stamp


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
        "matpredict_version": running_version(),
        # The checkout the code ran from (issue #139); absent when it cannot be read.
        **({"code": stamp} if (stamp := code_stamp()) else {}),
        "database": database_digest(db_root),
        "taxonomy_source": taxonomy_source,
        "genome": genome_doc,
        "proteins": Path(proteins).name if proteins else None,
        "parameters": parameters,
        "started": clock.started,
        "wall_seconds": clock.wall_seconds(),
    }
