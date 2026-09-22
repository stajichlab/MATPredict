"""taxonkit subprocess wrapper for lineage resolution, plus an NCBI-efetch-backed
ancestor-taxid lookup used for lineage-aware `taxonomic_scope` matching."""
from __future__ import annotations

import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


@dataclass(frozen=True)
class TaxonomyResult:
    """Result of resolving a taxid's lineage via taxonkit."""

    taxid: int
    lineage: str
    is_current: bool


def resolve_lineage(taxid: int, runner: Callable = subprocess.run) -> TaxonomyResult:
    """Resolve a taxid to its full lineage string via `taxonkit reformat`."""
    proc = runner(
        ["taxonkit", "reformat", "-I", "1", "-f", "k__{k};p__{p};c__{c};o__{o};f__{f};g__{g};s__{s}"],
        input=str(taxid),
        capture_output=True,
        text=True,
    )
    line = proc.stdout.strip().split("\n")[0] if proc.stdout.strip() else f"{taxid}\t"
    _, _, lineage = line.partition("\t")
    is_current = bool(lineage.strip())
    return TaxonomyResult(taxid=taxid, lineage=lineage, is_current=is_current)


def _default_transport(url: str, max_attempts: int = 4, backoff_seconds: float = 2.0) -> str:
    """GET a URL, retrying on rate-limit/server errors -- mirrors `db.cli._http_transport`
    but is defined locally to avoid a `taxonomy -> cli -> validate -> taxonomy` import
    cycle (`db.validate` already imports this module)."""
    import requests

    last_error: Exception | None = None
    for attempt in range(max_attempts):
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
        last_error = requests.HTTPError(f"{response.status_code} for {url}: {response.text[:200]}")
        if response.status_code in (429, 500, 502, 503, 504) and attempt < max_attempts - 1:
            time.sleep(backoff_seconds * (attempt + 1))
            continue
        break
    raise last_error


_default_ncbi_client = None


def _get_default_ncbi_client():
    """Lazily build a process-lifetime default `NcbiClient` for lineage lookups.

    Reuses `MATPredict.db.ncbi_client.NcbiClient` and `MATPredict.db.http_cache.CachedFetcher`
    -- the repo's existing NCBI client/on-disk-cache infrastructure -- rather than a
    separate mechanism. The `CachedFetcher` cache directory follows the same convention
    as `MatpredictConfig.from_env` (`.matpredict_cache` under the current working
    directory, or `$MATPREDICT_CACHE_DIR`), so lineage lookups persist across runs the
    same way accession/sequence lookups already do.
    """
    global _default_ncbi_client
    if _default_ncbi_client is None:
        import os

        from MATPredict.db.http_cache import CachedFetcher
        from MATPredict.db.ncbi_client import NcbiClient

        cache_dir = Path(os.environ.get("MATPREDICT_CACHE_DIR", Path.cwd() / ".matpredict_cache"))
        email = os.environ.get("MATPREDICT_NCBI_EMAIL", "jason.stajich@ucr.edu")
        api_key = os.environ.get("MATPREDICT_NCBI_API_KEY")
        fetcher = CachedFetcher(cache_dir=cache_dir, transport=_default_transport)
        _default_ncbi_client = NcbiClient(email=email, api_key=api_key, fetcher=fetcher)
    return _default_ncbi_client


def default_lineage_taxids(taxid: int) -> list[int]:
    """Fetch `taxid`'s ancestor lineage as a list of taxids using a lazily-built default
    `NcbiClient` (see `_get_default_ncbi_client`). This is the default
    `lineage_taxids_resolver` for `detect.family_registry.route`."""
    return _get_default_ncbi_client().fetch_taxonomy_lineage(taxid)


def default_lineage_phylum_name(taxid: int) -> str | None:
    """Return the scientific name of `taxid`'s phylum-rank ancestor, or None.

    This is the default `phylum_name_resolver` for
    `detect.family_registry.route`'s phylum fallback. It reuses the same lazily
    built default `NcbiClient` -- and therefore the same on-disk, URL-keyed
    response cache -- as `default_lineage_taxids`, and both read the identical
    `efetch db=taxonomy&id=<taxid>` document, so asking for the phylum after
    asking for the lineage of the same taxid is a cache hit, not a second
    network request.
    """
    return _get_default_ncbi_client().fetch_taxonomy_phylum(taxid)


def default_genetic_code(taxid: int) -> int | None:
    """The NCBI translation table for `taxid`, or None when unknown.

    Same lazily-built client and the same on-disk URL-keyed cache as
    `default_lineage_taxids`, reading the identical efetch document, so a run
    that already routed by taxid pays nothing extra for this.
    """
    return _get_default_ncbi_client().fetch_taxonomy_genetic_code(taxid)
