"""Runtime configuration for MATPredict's db subcommands."""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MatpredictConfig:
    """Resolved configuration for a matpredict invocation."""

    db_root: Path
    ncbi_email: str
    ncbi_api_key: str | None
    cache_dir: Path

    @classmethod
    def from_env(cls, repo_root: Path) -> "MatpredictConfig":
        """Build config from the repo layout plus environment overrides."""
        db_root = Path(os.environ.get("MATPREDICT_DB_ROOT", repo_root / "db"))
        ncbi_email = os.environ.get("MATPREDICT_NCBI_EMAIL", "jason.stajich@ucr.edu")
        ncbi_api_key = os.environ.get("MATPREDICT_NCBI_API_KEY")
        cache_dir = Path(os.environ.get("MATPREDICT_CACHE_DIR", repo_root / ".matpredict_cache"))
        return cls(db_root=db_root, ncbi_email=ncbi_email, ncbi_api_key=ncbi_api_key, cache_dir=cache_dir)
