"""scripts/annotation_by_tier.py joins verdict rows to tiers by the whole-list locus index."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import yaml

SPEC = importlib.util.spec_from_file_location(
    "annotation_by_tier", Path(__file__).resolve().parents[1] / "scripts" / "annotation_by_tier.py")
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def _run(tmp_path):
    d = tmp_path / "runs" / "S1"
    d.mkdir(parents=True)
    (d / "detection_report.yaml").write_text(yaml.safe_dump({"detected": [
        {"family": "Basidiomycota:HD", "confidence": "high", "idiomorph": "undetermined"},
        {"family": "Basidiomycota:PR", "confidence": "medium", "idiomorph": "undetermined"}]}))
    v = tmp_path / "v.tsv"
    v.write_text("sample\tlocus\tmp_gene\tverdict\n"
                 "S1\tBasidiomycota_HD_0\tHD1\tannotated\n"
                 "S1\tBasidiomycota_PR_1\tSTE3\tnovel\n"
                 "S1\tBasidiomycota_PR_0\tX\tpartial\n")  # PR is index 1 overall, so _0 does not match
    return tmp_path / "runs", v


def test_index_is_the_position_in_the_whole_detected_list(tmp_path):
    runs, v = _run(tmp_path)
    rows, counts = mod.join(runs, v)
    assert [r["tier"] for r in rows] == ["high", "medium", "unmatched"]
    assert counts["medium"]["novel"] == 1 and counts["unmatched"]["partial"] == 1


def test_main_writes_both_tables(tmp_path):
    runs, v = _run(tmp_path)
    assert mod.main([str(runs), str(v), str(tmp_path / "out")]) == 0
    assert (tmp_path / "out_by_tier.tsv").read_text().splitlines()[1].startswith("high\t1\t1\t0\t0")
    assert (tmp_path / "out_rows.tsv").exists()
