from __future__ import annotations

from MATPredict.__main__ import main


def test_help_exits_zero(capsys):
    exit_code = main(["curate-db", "--help"])
    assert exit_code == 0
    captured = capsys.readouterr()
    assert "curate-db" in captured.out or "usage" in captured.out
