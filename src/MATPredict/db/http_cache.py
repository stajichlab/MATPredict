"""On-disk HTTP response cache, keyed by URL hash, for NCBI/UniProt clients."""
from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


@dataclass
class CachedFetcher:
    """Caches transport(url) results as files under cache_dir, keyed by url hash."""

    cache_dir: Path
    transport: Callable[[str], str]

    def _cache_path(self, url: str) -> Path:
        digest = hashlib.sha256(url.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{digest}.txt"

    def get(self, url: str) -> str:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        path = self._cache_path(url)
        # An empty file is a miss. Before writes were atomic, a reader racing
        # a writer could see a truncated file and return '' as the answer.
        if path.exists() and path.stat().st_size > 0:
            return path.read_text()
        body = self.transport(url)
        # Write to a unique temporary name, then rename into place: a rename
        # is atomic within a filesystem, so a concurrent reader (many detect
        # processes share one cache) sees nothing or the whole response.
        tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
        tmp.write_text(body)
        os.replace(tmp, path)
        return body
