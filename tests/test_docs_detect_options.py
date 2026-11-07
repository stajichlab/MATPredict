"""docs/detect-options.md must name every policy option, preset and environment variable the code defines (it is the one page that lists them)."""
import re
from pathlib import Path

from MATPredict.__main__ import build_parser
from MATPredict.detect.cli import PRESETS

DOC = (Path(__file__).resolve().parent.parent / "docs" / "detect-options.md").read_text()
SRC = Path(__file__).resolve().parent.parent / "src" / "MATPredict"


def _detect_parser():
    sub = build_parser()._subparsers._group_actions[0].choices
    return sub["detect"]


def test_every_detect_policy_option_is_documented():
    names = [a.option_strings[0] for a in _detect_parser()._actions if a.option_strings and a.dest.endswith("_policy")]
    assert names, "no policy options found"
    assert [n for n in names if f"`{n}`" not in DOC] == []


def test_every_preset_is_documented():
    assert [p for p in PRESETS if f"`{p}`" not in DOC and f"`--preset {p}`" not in DOC] == []


def test_every_environment_variable_is_documented():
    used = set()
    for path in SRC.rglob("*.py"):
        used.update(re.findall(r"MATPREDICT_[A-Z0-9_]+", path.read_text()))
    assert sorted(v for v in used if v not in DOC) == []
