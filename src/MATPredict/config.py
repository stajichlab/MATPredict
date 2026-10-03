"""Runtime configuration for MATPredict's db subcommands."""
from __future__ import annotations

import os
import tomllib
from dataclasses import dataclass
from pathlib import Path

#: Config file read when $MATPREDICT_CONFIG is not set.
DEFAULT_CONFIG_PATH = Path("~/.config/matpredict/config.toml")


def config_path() -> Path:
    """$MATPREDICT_CONFIG, else $XDG_CONFIG_HOME/matpredict/config.toml, else
    ~/.config/matpredict/config.toml."""
    explicit = os.environ.get("MATPREDICT_CONFIG")
    if explicit:
        return Path(explicit).expanduser()
    xdg = os.environ.get("XDG_CONFIG_HOME")
    if xdg:
        return Path(xdg) / "matpredict" / "config.toml"
    return DEFAULT_CONFIG_PATH.expanduser()


def _read_config_file() -> dict:
    path = config_path()
    if not path.is_file():
        return {}
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def ncbi_identity() -> tuple[str | None, str | None]:
    """(email, api_key) sent to NCBI E-utilities.

    Each value comes from the environment ($MATPREDICT_NCBI_EMAIL,
    $MATPREDICT_NCBI_API_KEY) or, when that is unset, from the `[ncbi]` table of
    the config file (`email`, `api_key`). There is no default: with no e-mail,
    requests carry none, and `NcbiClient` logs one warning per process at the
    first request (not at start-up, so `--help` and offline runs stay quiet).
    """
    ncbi = (_read_config_file().get("ncbi") or {})
    email = os.environ.get("MATPREDICT_NCBI_EMAIL") or ncbi.get("email") or None
    api_key = os.environ.get("MATPREDICT_NCBI_API_KEY") or ncbi.get("api_key") or None
    return email, api_key


@dataclass(frozen=True)
class MatpredictConfig:
    """Resolved configuration for a matpredict invocation."""

    db_root: Path
    ncbi_email: str | None
    ncbi_api_key: str | None
    cache_dir: Path

    @classmethod
    def from_env(cls, repo_root: Path) -> "MatpredictConfig":
        """Build config from the repo layout plus environment overrides."""
        db_root = Path(os.environ.get("MATPREDICT_DB_ROOT", repo_root / "db"))
        ncbi_email, ncbi_api_key = ncbi_identity()
        cache_dir = Path(os.environ.get("MATPREDICT_CACHE_DIR", repo_root / ".matpredict_cache"))
        return cls(db_root=db_root, ncbi_email=ncbi_email, ncbi_api_key=ncbi_api_key, cache_dir=cache_dir)
