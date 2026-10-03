"""On-disk HTTP response cache, keyed by URL hash, for NCBI/UniProt clients."""
from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

#: Query parameters that identify the caller, not the request. They are left out
#: of the cache key, so the cache is shared whatever e-mail or API key is set.
IDENTITY_PARAMS = ("email", "api_key", "tool")

#: Until 2026-10-03 the cache key was the full URL, and every NCBI URL ended with
#: this fixed default e-mail. A miss on the new key retries the old key once and
#: copies the entry forward, so caches built before then stay valid. The
#: address is used only to find old cache files; it is never sent to NCBI.
_LEGACY_EMAIL_SUFFIX = "&email=jason.stajich@ucr.edu"


def cache_key(url: str) -> str:
    """`url` without the caller-identity query parameters (`IDENTITY_PARAMS`).

    Plain string edit, no re-encoding: every other byte of the URL is kept, so the
    key of a URL with no identity parameters is the URL itself (as before 2026-10-03).
    """
    base, sep, query = url.partition("?")
    if not sep:
        return url
    kept = [item for item in query.split("&") if item.split("=", 1)[0] not in IDENTITY_PARAMS]
    return f"{base}?{'&'.join(kept)}" if kept else base


@dataclass
class CachedFetcher:
    """Caches transport(url) results as files under cache_dir, keyed by url hash."""

    cache_dir: Path
    transport: Callable[[str], str]

    def _raw_path(self, key: str) -> Path:
        digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{digest}.txt"

    def _cache_path(self, url: str) -> Path:
        return self._raw_path(cache_key(url))

    def _legacy_path(self, url: str) -> Path | None:
        """The pre-2026-10-03 cache file of an NCBI URL (see `_LEGACY_EMAIL_SUFFIX`)."""
        key = cache_key(url)
        if "eutils.ncbi.nlm.nih.gov" not in key:
            return None
        return self._raw_path(key + _LEGACY_EMAIL_SUFFIX)

    def get(self, url: str) -> str:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        path = self._cache_path(url)
        # An empty file is a miss. Before writes were atomic, a reader racing
        # a writer could see a truncated file and return '' as the answer.
        if path.exists() and path.stat().st_size > 0:
            return path.read_text()
        legacy = self._legacy_path(url)
        if legacy is not None and legacy.exists() and legacy.stat().st_size > 0:
            body = legacy.read_text()
            self._write(path, body)
            return body
        body = self.transport(url)
        # Write to a unique temporary name, then rename into place: a rename
        # is atomic within a filesystem, so a concurrent reader (many detect
        # processes share one cache) sees nothing or the whole response.
        self._write(path, body)
        return body

    @staticmethod
    def _write(path: Path, body: str) -> None:
        tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
        tmp.write_text(body)
        os.replace(tmp, path)
