"""NCBI identity: environment first, then the config file, and no default."""
import logging

import pytest

import MATPredict.config as config
from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import NcbiClient


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch, tmp_path):
    for v in ("MATPREDICT_NCBI_EMAIL", "MATPREDICT_NCBI_API_KEY", "MATPREDICT_CONFIG", "XDG_CONFIG_HOME"):
        monkeypatch.delenv(v, raising=False)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setattr(config, "_warned_no_email", False)


def test_no_default_email(caplog):
    with caplog.at_level(logging.WARNING):
        assert config.ncbi_identity() == (None, None)
    assert "no NCBI e-mail configured" in caplog.text


def test_the_config_file_supplies_email_and_key(monkeypatch, tmp_path):
    f = tmp_path / "c.toml"
    f.write_text('[ncbi]\nemail = "file@example.org"\napi_key = "K1"\n')
    monkeypatch.setenv("MATPREDICT_CONFIG", str(f))
    assert config.ncbi_identity() == ("file@example.org", "K1")


def test_the_default_config_path_is_under_home(tmp_path):
    p = tmp_path / "home" / ".config" / "matpredict" / "config.toml"
    p.parent.mkdir(parents=True)
    p.write_text('[ncbi]\nemail = "home@example.org"\n')
    assert config.ncbi_identity() == ("home@example.org", None)


def test_xdg_config_home_is_honoured(monkeypatch, tmp_path):
    p = tmp_path / "xdg" / "matpredict" / "config.toml"
    p.parent.mkdir(parents=True)
    p.write_text('[ncbi]\nemail = "xdg@example.org"\n')
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))
    assert config.ncbi_identity()[0] == "xdg@example.org"


def test_the_environment_wins_over_the_file(monkeypatch, tmp_path):
    f = tmp_path / "c.toml"
    f.write_text('[ncbi]\nemail = "file@example.org"\napi_key = "K1"\n')
    monkeypatch.setenv("MATPREDICT_CONFIG", str(f))
    monkeypatch.setenv("MATPREDICT_NCBI_EMAIL", "env@example.org")
    assert config.ncbi_identity() == ("env@example.org", "K1")


def test_a_request_without_email_carries_none(tmp_path):
    seen = []
    client = NcbiClient(email=None, api_key=None,
                        fetcher=CachedFetcher(cache_dir=tmp_path, transport=lambda u: seen.append(u) or "<x/>"))
    client.fetcher.get(client._url("efetch.fcgi", "db=taxonomy&id=4837&retmode=xml"))
    assert "email=" not in seen[0] and "tool=MATPredict" in seen[0]


def test_a_configured_email_is_sent():
    client = NcbiClient(email="me@example.org", api_key=None, fetcher=None)
    assert "&email=me@example.org" in client._url("efetch.fcgi", "db=taxonomy&id=1")
