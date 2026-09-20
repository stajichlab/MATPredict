"""Shared pytest fixtures.

`real_db_root` exists so tests that measure something against the LIVE curated
database do not each hardcode `Path("db")`, which silently depends on pytest
being invoked from the repository root and gives a confusing "no such file"
failure when it is not. The path is resolved from this file's own location
instead, so it is correct from any working directory.
"""
from __future__ import annotations

from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def real_db_root() -> Path:
    """The repository's live `db/` directory, skipping the test if it is absent."""
    db_root = _REPO_ROOT / "db"
    if not db_root.is_dir():
        pytest.skip(f"the live curated database is not present at {db_root}")
    return db_root
