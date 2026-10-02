"""scripts/submit_binned_panel.py: genome-size bins (curator ruling 2026-10-01).

Loaded by path because `scripts/` is not an importable package.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "submit_binned_panel.py"


def _load():
    spec = importlib.util.spec_from_file_location("submit_binned_panel", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_large_genomes_go_to_their_own_bin():
    mod = _load()
    samples = [("a", "1"), ("b", "2"), ("big", "3")]
    sizes = {"a": 40_000_000, "b": 90_000_000, "big": 1_710_000_000}
    chunks, _, no_size = mod.plan(samples, sizes, large_bp=500_000_000, small_cpus=16,
                                  small_seconds=330, target_hours=1.25, large_per_job=64)
    assert [(c["name"], [r[0] for r in c["rows"]]) for c in chunks] == [
        ("small_000", ["a", "b"]), ("large_000", ["big"])]
    assert no_size == 0


def test_small_chunk_size_follows_target_runtime():
    mod = _load()
    # 1.25 h * 3600 * 16 / 330 s = 218.18 -> 218 genomes per chunk
    samples = [(f"g{i}", "1") for i in range(500)]
    chunks, per_chunk, _ = mod.plan(samples, {}, large_bp=500_000_000, small_cpus=16,
                                    small_seconds=330, target_hours=1.25, large_per_job=64)
    assert per_chunk == 218
    assert [len(c["rows"]) for c in chunks] == [218, 218, 64]


def test_missing_size_counts_and_goes_small():
    mod = _load()
    chunks, _, no_size = mod.plan([("x", "1")], {}, large_bp=500_000_000, small_cpus=16,
                                  small_seconds=330, target_hours=1.25, large_per_job=64)
    assert no_size == 1 and chunks[0]["bin"] == "small"


def test_clade_panel_reads_genome_timeout():
    text = (SCRIPT.parent / "run_clade_panel.slurm").read_text()
    assert 'GENOME_TIMEOUT="${GENOME_TIMEOUT:-3600}"' in text
    assert 'timeout "$GENOME_TIMEOUT"' in text
    assert "timeout 3600" not in text
