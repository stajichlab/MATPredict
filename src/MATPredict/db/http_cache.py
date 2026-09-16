"""On-disk HTTP response cache, keyed by URL hash, for NCBI/UniProt clients."""
from __future__ import annotations

import hashlib
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
        if path.exists():
            return path.read_text()
        body = self.transport(url)
        path.write_text(body)
        return body
