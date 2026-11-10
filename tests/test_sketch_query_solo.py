"""datasets/reads/sketch_store/sketch_query.py solo_hashes: a 'specific' hash that keeps turning up without the rest of its target is dropped."""
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("sketch_query", Path(__file__).resolve().parents[1] / "datasets" / "reads" / "sketch_store" / "sketch_query.py")
sq = importlib.util.module_from_spec(spec); spec.loader.exec_module(sq)


def test_a_hash_seen_alone_in_many_samples_is_dropped():
    # bafA-like: hashes 1, 2, 3; 40 samples carry the gene (all three), 300 have hash 1 alone, 660 have none
    samples = [{1, 2, 3}] * 40 + [{1}] * 300 + [set()] * 660
    assert sq.solo_hashes({1, 2, 3}, samples) == {1: 300}


def test_rare_solo_hits_and_single_hash_targets_are_kept():
    samples = [{1, 2, 3}] * 40 + [{2}] * 5 + [set()] * 955        # 5 solo hits: under max(10, 1%)
    assert sq.solo_hashes({1, 2, 3}, samples) == {}
    assert sq.solo_hashes({1}, [{1}] * 500) == {}                  # one specific hash: nothing to compare it with
    assert sq.solo_hashes({1, 2, 3}, [{1}] * 500, max_solo=0) == {}  # switched off


def test_every_hash_of_a_target_is_never_dropped():
    samples = [{1}] * 200 + [{2}] * 200
    assert sq.solo_hashes({1, 2}, samples) == {}
