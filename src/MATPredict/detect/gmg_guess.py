"""Guess a genome's phylum with GmG ("Guess My Genome") so `matpredict detect --phylum auto` can narrow the search.

GmG (https://github.com/stajichlab/GmG) compares a sketch of the genome with a reference index of fungal genomes and calls
the taxonomy from phylum to species. This module only (1) finds the index, (2) runs GmG, (3) turns its record into a
`GmgGuess`, and (4) decides what MATPredict does with it (`decide`). Routing itself stays in `family_registry.route`.

Design: docs/notes/2026-11-02_gmg-phylum-routing.md. In short, a wrong phylum is worse than none (the right family is never
searched), so a call that GmG itself marks weak, or an input it says is not a fungus, is not used unless the operator asks
(`--phylum auto-weak`), and the run then says plainly that it was not searched and why.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

from MATPredict import logger

#: Where GmG publishes its indexes for the web (CORS-open, also fine for the CLI).
GMG_INDEX_SITE = "https://gmg.fungalgenomes.org/index/"
#: Release assets (and SHA256SUMS) of one index release: `<base>/index-<version>/<file>`.
GMG_RELEASE_BASE = "https://github.com/stajichlab/GmG/releases/download/"
#: The index MATPredict is validated against (decision 2 of the design note: pin, and warn when a newer one exists).
DEFAULT_GMG_INDEX = "bfd-reps-2026.10.5-s20000.gmgi"
#: `--phylum` values that mean "ask GmG". `auto-weak` also accepts a call GmG marked `weak_phylum`.
AUTO = "auto"
AUTO_WEAK = "auto-weak"
AUTO_CHOICES = (AUTO, AUTO_WEAK)

#: File extensions `gmg run` finds as sequence files (plain or .gz/.zst/.bz2/.xz). A genome with another extension (`.fsa`, the GenBank
#: submission name, and `.txt` are common) is handed to GmG through a link named `.fna`: the content is the same.
GMG_SEQUENCE_EXTENSIONS = (".fa", ".fna", ".fasta", ".fas", ".ffn", ".fq", ".fastq")
_COMPRESSION_SUFFIXES = (".gz", ".zst", ".bz2", ".xz")
_HTTP_TIMEOUT = 30
_MANIFEST_TTL_SECONDS = 24 * 3600
_VERSION_RE = re.compile(r"(\d{4}\.\d+(?:\.\d+)?)")


class GmgUnavailable(RuntimeError):
    """GmG (Python bindings or binary) or its index cannot be used. A usage error: it says what to install."""


@dataclass
class GmgGuess:
    """What GmG said about one genome, in the terms `route()` needs."""

    status: str | None = None             # GmG record status: ok | low_confidence | no_match | error
    usable: bool = False                  # GmG's own routing_hint.usable
    reason: str | None = None             # why not usable: no_phylum_call | not_fungal | weak_phylum | error
    phylum: str | None = None
    phylum_taxid: int | None = None
    deepest_rank: str | None = None
    #: NCBI taxid of the family call, else the order call, else None. Never genus or species: with a new genus the
    #: genus taxid is the nearest relative's, not the genome's (measured on the regression panel, design note).
    routing_taxid: int | None = None
    routing_taxid_rank: str | None = None
    flags: list[str] = field(default_factory=list)
    gmg_version: str | None = None
    index: str | None = None              # index file name
    index_sha256: str | None = None
    index_spec: str | None = None         # what the operator asked for (a path, a file name, `latest` or a URL)

    def to_dict(self) -> dict:
        return asdict(self)


# --------------------------------------------------------------------------------------------- the GmG record -> a guess

def guess_from_record(rec: dict) -> GmgGuess:
    """Read the `routing_hint` of a GmG record (schema 9). A record of an older GmG has none, so the same rule is applied
    here: usable = a phylum is called, the input is not `not_fungal`, and the call is not `weak_phylum`."""
    flags = list(rec.get("qc_flags") or [])
    calls = ((rec.get("guess") or {}).get("calls")) or {}

    def call(rank: str) -> dict:
        return calls.get(rank) or {}

    hint = rec.get("routing_hint")
    if hint is None:
        phylum = call("phylum").get("taxon")
        reason = ("error" if rec.get("status") == "error" else
                  "no_phylum_call" if not phylum else
                  "not_fungal" if "not_fungal" in flags else
                  "weak_phylum" if "weak_phylum" in flags else None)
        hint = {"usable": reason is None, "reason": reason, "phylum": phylum,
                "phylum_taxid": call("phylum").get("taxid")}
    routing_taxid = rank = None
    for r in ("family", "order"):
        tid = call(r).get("taxid")
        if call(r).get("taxon") and tid:
            routing_taxid, rank = int(tid), r
            break
    deepest = hint.get("deepest_rank")
    if deepest is None:
        called = [r for r in ("phylum", "class", "order", "family", "genus", "species") if call(r).get("taxon")]
        deepest = called[-1] if called else None
    return GmgGuess(
        status=rec.get("status"), usable=bool(hint.get("usable")), reason=hint.get("reason"),
        phylum=hint.get("phylum"), phylum_taxid=hint.get("phylum_taxid"), deepest_rank=deepest,
        routing_taxid=routing_taxid, routing_taxid_rank=rank, flags=flags,
    )


# ------------------------------------------------------------------------------------------------------ the decision

@dataclass(frozen=True)
class GmgDecision:
    """What `--phylum auto` does with a guess.

    `phylum` is the curated phylum to route to, or None for "do not search". `note` is one sentence for the report either
    way. `weak` is True when a `weak_phylum` call was accepted (`auto-weak`); the report then says so."""

    phylum: str | None
    taxid: int | None
    note: str
    weak: bool = False


def decide(guess: GmgGuess, available: list[str], accept_weak: bool = False) -> GmgDecision:
    """Turn a guess into a routing decision. `available` is `available_phyla(db_root)`, the phyla with curated families.

    * usable and in `available` -> route to it (with the order/family taxid when GmG gave one);
    * `weak_phylum` and `accept_weak` and in `available` -> route to it, marked weak;
    * a fungal phylum with no curated family -> not searched (the same outcome a taxid in that phylum gives today);
    * everything else (no call, `not_fungal`, a weak call without `accept_weak`, an error) -> not searched, with the reason.
    """
    ph = guess.phylum
    if guess.reason == "error":
        return GmgDecision(None, None, "GmG could not read the genome; not searched")
    if not ph or guess.reason == "no_phylum_call":
        return GmgDecision(None, None, "GmG made no phylum call for this genome; not searched")
    if guess.reason == "not_fungal":
        return GmgDecision(None, None, f"GmG says this is probably not a fungus (nearest phylum {ph}); not searched")
    weak = guess.reason == "weak_phylum"
    if weak and not accept_weak:
        return GmgDecision(None, None, f"GmG's phylum call {ph} is weak (best-hit ANI below 0.85); not searched "
                                       f"(use --phylum {AUTO_WEAK} to search it anyway)")
    if guess.reason not in (None, "weak_phylum"):  # a reason added by a newer GmG: be cautious
        return GmgDecision(None, None, f"GmG marks its call {ph} as not usable ({guess.reason}); not searched")
    if ph not in available:
        return GmgDecision(None, None, f"GmG calls phylum {ph}, which has no curated MAT locus family; not searched")
    note = f"GmG called phylum {ph}" + (" (weak call, accepted by --phylum auto-weak)" if weak else "")
    return GmgDecision(ph, guess.routing_taxid, note, weak=weak)


# ------------------------------------------------------------------------------------------------------- the index

def gmg_cache_dir() -> Path:
    """Where the GmG index is kept: `$MATPREDICT_GMG_CACHE`, else `$MATPREDICT_CACHE_DIR/gmg`, else
    `$XDG_CACHE_HOME/matpredict/gmg`, else `~/.cache/matpredict/gmg` (decision 1 of the design note)."""
    if os.environ.get("MATPREDICT_GMG_CACHE"):
        return Path(os.environ["MATPREDICT_GMG_CACHE"]).expanduser()
    if os.environ.get("MATPREDICT_CACHE_DIR"):
        return Path(os.environ["MATPREDICT_CACHE_DIR"]).expanduser() / "gmg"
    base = Path(os.environ["XDG_CACHE_HOME"]) if os.environ.get("XDG_CACHE_HOME") else Path("~/.cache").expanduser()
    return base / "matpredict" / "gmg"


def _get(url: str, **kw):
    import requests
    return requests.get(url, timeout=_HTTP_TIMEOUT, **kw)


def _manifest(cache: Path) -> list[dict] | None:
    """The site's `manifest.json`, cached for a day; None when it cannot be had (never raises)."""
    path = cache / "manifest.json"
    try:
        if path.is_file() and time.time() - path.stat().st_mtime < _MANIFEST_TTL_SECONDS:
            return json.loads(path.read_text())
        r = _get(GMG_INDEX_SITE + "manifest.json")
        r.raise_for_status()
        cache.mkdir(parents=True, exist_ok=True)
        path.write_text(r.text)
        return r.json()
    except Exception as exc:  # noqa: BLE001 - the manifest is an aid (latest, newer-index warning), not a requirement
        logger.debug("GmG index manifest not available: %s: %s", type(exc).__name__, exc)
        try:
            return json.loads(path.read_text()) if path.is_file() else None
        except Exception:  # noqa: BLE001
            return None


def _version(name: str) -> tuple[int, ...]:
    m = _VERSION_RE.search(name)
    return tuple(int(x) for x in m.group(1).split(".")) if m else ()


def newer_index(name: str, manifest: list[dict] | None) -> str | None:
    """The newest representatives index in `manifest` when it is newer than `name`, else None."""
    if not manifest or not _version(name):
        return None
    reps = [m["file"] for m in manifest if "reps" in m.get("file", "") and _version(m["file"])]
    best = max(reps, key=_version, default=None)
    return best if best and _version(best) > _version(name) else None


def latest_index(manifest: list[dict] | None) -> str:
    reps = [m["file"] for m in (manifest or []) if "reps" in m.get("file", "") and _version(m["file"])]
    if not reps:
        raise GmgUnavailable(f"could not find the newest GmG index at {GMG_INDEX_SITE}manifest.json; "
                             "pass --gmg-index with a file name or a path")
    return max(reps, key=_version)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _release_sha256(name: str) -> str | None:
    """The sha256 GmG published for `name` in its release's SHA256SUMS, or None when that cannot be read."""
    m = _VERSION_RE.search(name)
    if not m:
        return None
    try:
        r = _get(f"{GMG_RELEASE_BASE}index-{m.group(1)}/SHA256SUMS")
        r.raise_for_status()
        for line in r.text.splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[1].lstrip("*") == name:
                return parts[0]
    except Exception as exc:  # noqa: BLE001
        logger.debug("SHA256SUMS for %s not available: %s: %s", name, type(exc).__name__, exc)
    return None


def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=dest.parent, suffix=".part")
    try:
        with os.fdopen(fd, "wb") as out, _get(url, stream=True) as r:
            r.raise_for_status()
            for chunk in r.iter_content(1 << 20):
                out.write(chunk)
        os.replace(tmp, dest)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def resolve_index(spec: str | None, *, cache: Path | None = None) -> tuple[Path, str]:
    """Find (downloading once if needed) the GmG index, returning (path, the spec used).

    `spec`: an existing file path; `latest` (the newest representatives index on the GmG site); a release file name such as
    `bfd-reps-2026.10.5-s20000.gmgi`; or an http(s) URL. None means `$MATPREDICT_GMG_INDEX`, else the pinned default.
    A downloaded file is checked against the sha256 in GmG's release when that can be read, and kept in the cache.
    When the index is a pinned file name and a newer one exists, one warning says so; the answer never changes silently."""
    spec = spec or os.environ.get("MATPREDICT_GMG_INDEX") or DEFAULT_GMG_INDEX
    if spec.startswith(("http://", "https://")):
        url, name = spec, Path(spec.split("?")[0]).name
    elif Path(spec).expanduser().is_file():
        return Path(spec).expanduser().resolve(), spec
    elif "/" in spec or spec.endswith(os.sep):
        raise GmgUnavailable(f"GmG index not found: {spec}")
    else:
        cache = cache or gmg_cache_dir()
        manifest = _manifest(cache)
        name = latest_index(manifest) if spec == "latest" else spec
        url = GMG_INDEX_SITE + name
        if spec != "latest":
            newer = newer_index(name, manifest)
            if newer:
                logger.warning("a newer GmG index exists (%s); this run uses the pinned %s. Pass --gmg-index %s "
                               "or set MATPREDICT_GMG_INDEX to change it", newer, name, newer)
    cache = cache or gmg_cache_dir()
    dest = cache / name
    if not dest.is_file():
        logger.info("downloading the GmG index %s to %s (one time)", url, cache)
        try:
            _download(url, dest)
        except Exception as exc:  # noqa: BLE001
            raise GmgUnavailable(f"could not download the GmG index {url}: {type(exc).__name__}: {exc}. "
                                 f"Download it yourself and pass --gmg-index PATH") from exc
        want = _release_sha256(name)
        if want and _sha256(dest) != want:
            dest.unlink(missing_ok=True)
            raise GmgUnavailable(f"the downloaded GmG index {name} does not match its published sha256; "
                                 "removed it, try again")
    return dest, spec


# --------------------------------------------------------------------------------------------------------- run GmG

def _as_gmg_input(genome: str | Path, tmp: Path) -> Path:
    """The genome as `gmg run` accepts it: unchanged when its extension is one GmG recognises, else a link named `genome.fna` in `tmp`
    (GmG finds sequence files by extension and refuses `.fsa`; the Python bindings read any path, so only the binary needs this)."""
    path = Path(genome)
    name = path.name.lower()
    for suffix in _COMPRESSION_SUFFIXES:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    if name.endswith(GMG_SEQUENCE_EXTENSIONS):
        return path
    link = tmp / ("genome.fna" + next((x for x in _COMPRESSION_SUFFIXES if path.name.lower().endswith(x)), ""))
    link.symlink_to(path.resolve())
    return link


def classify_genome(genome: str | Path, index_path: Path) -> tuple[dict, str | None]:
    """Run GmG on one genome and return (its record, the GmG version). Uses GmG's Python bindings when they are installed,
    else the `gmg` binary on PATH (`gmg run`). Raises `GmgUnavailable` when neither is there.

    NOTE: the name `gmg` on PyPI belongs to an unrelated project, so GmG's wheel is not installed with `pip install gmg`
    (its distribution is `guess-my-genome`, not yet on PyPI), and an `import gmg` that has no `Index` is that other package and is ignored."""
    try:
        import gmg  # type: ignore[import-not-found]
    except ImportError:
        gmg = None
    if gmg is not None and not hasattr(gmg, "Index"):
        gmg = None
    if gmg is not None:
        rec = gmg.Index(str(index_path)).classify(str(genome))
        return rec, getattr(gmg, "__version__", None)
    exe = shutil.which("gmg")
    if not exe:
        raise GmgUnavailable("--phylum auto needs GmG: put the `gmg` binary on PATH (release binaries for linux-64 and osx-arm64: "
                             "https://github.com/stajichlab/GmG/releases) or install GmG's Python bindings (distribution "
                             "`guess-my-genome`, built from that repository; not `pip install gmg`, which is a different project). Or give --phylum or --taxid")
    with tempfile.TemporaryDirectory(prefix="matpredict_gmg_") as tmp:
        proc = subprocess.run([exe, "-q", "run", str(_as_gmg_input(genome, Path(tmp))), "--index", str(index_path), "-o", tmp + "/out", "--no-report"],
                              capture_output=True, text=True)
        out = Path(tmp) / "out" / "results.jsonl"
        if proc.returncode != 0 or not out.is_file():
            raise GmgUnavailable(f"`gmg run` failed (exit {proc.returncode}): {(proc.stderr or '').strip()[:300]}")
        rec = json.loads(out.read_text().splitlines()[0])
    ver = subprocess.run([exe, "--version"], capture_output=True, text=True).stdout.strip() or None
    return rec, ver


def guess_phylum(genome: str | Path, index_spec: str | None = None) -> GmgGuess:
    """Resolve the index, run GmG on `genome`, and return the `GmgGuess` (with the GmG and index identity filled in)."""
    index_path, spec = resolve_index(index_spec)
    rec, version = classify_genome(genome, index_path)
    g = guess_from_record(rec)
    g.gmg_version, g.index, g.index_spec = version, index_path.name, spec
    g.index_sha256 = _sha256(index_path)
    return g
