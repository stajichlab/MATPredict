from __future__ import annotations

from types import SimpleNamespace

import pytest

from MATPredict.db.taxonomy import TaxonkitError, resolve_lineage


def _fake_runner_current(cmd, **kwargs):
    return SimpleNamespace(returncode=0, stdout="4837\tk__Fungi;p__Mucoromycota;...;s__Phycomyces_blakesleeanus\n")


def _fake_runner_merged(cmd, **kwargs):
    # taxonkit reformat prints an empty lineage for a taxid it can't resolve directly
    return SimpleNamespace(returncode=0, stdout="4837\t\n")


def test_resolve_current_taxid():
    result = resolve_lineage(4837, runner=_fake_runner_current)
    assert result.taxid == 4837
    assert "Mucoromycota" in result.lineage
    assert result.is_current is True


def test_resolve_merged_or_unknown_taxid():
    result = resolve_lineage(4837, runner=_fake_runner_merged)
    assert result.is_current is False
    assert result.lineage == ""


def test_the_transport_outlasts_a_burst_of_rate_limiting(monkeypatch):
    """128 concurrent runs share NCBI's keyless 3 req/s. Four attempts over
    ~12 s gave up inside a burst of 429s, and each give-up silently widened a
    run's routing. The default now rides out a longer burst."""
    import requests
    from MATPredict.db import taxonomy

    codes = iter([429] * 6 + [200])

    class R:
        def __init__(self, code):
            self.status_code, self.text = code, "ok" if code == 200 else "busy"

    monkeypatch.setattr(requests, "get", lambda url, **kw: R(next(codes)))
    monkeypatch.setattr(taxonomy.time, "sleep", lambda s: None)
    assert taxonomy._default_transport("https://example.org/t") == "ok"


def test_the_taxid_is_sent_newline_terminated():
    seen = {}

    def runner(cmd, **kwargs):
        seen.update(kwargs)
        return SimpleNamespace(returncode=0, stdout="5158\tk__;p__Ascomycota\n", stderr="")

    resolve_lineage(5158, runner=runner)
    assert seen["input"] == "5158\n"


def test_a_taxonkit_failure_is_an_error_not_an_out_of_date_taxid():
    def failing(cmd, **kwargs):
        return SimpleNamespace(returncode=1, stdout="", stderr="[ERRO] taxdump not found")

    with pytest.raises(TaxonkitError, match="taxdump not found"):
        resolve_lineage(5158, runner=failing)

    def silent(cmd, **kwargs):
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    with pytest.raises(TaxonkitError):
        resolve_lineage(5158, runner=silent)

    def missing(cmd, **kwargs):
        raise FileNotFoundError("taxonkit")

    with pytest.raises(TaxonkitError, match="not on PATH"):
        resolve_lineage(5158, runner=missing)
