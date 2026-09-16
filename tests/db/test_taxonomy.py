from __future__ import annotations

from types import SimpleNamespace

from MATPredict.db.taxonomy import resolve_lineage


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
