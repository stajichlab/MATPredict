from __future__ import annotations

from MATPredict.__main__ import main


def test_help_exits_zero(capsys):
    exit_code = main(["curate-db", "--help"])
    assert exit_code == 0
    captured = capsys.readouterr()
    # Verify action-specific content from curate-db subparser (not just top-level help)
    # Actions are registered in db/cli.py: propose, validate, accept, reject, build-gff, build-duckdb, release
    assert "propose" in captured.out or "validate" in captured.out
