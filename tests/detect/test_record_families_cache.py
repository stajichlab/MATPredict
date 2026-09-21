"""`load_record_families` parses every curated metadata.yaml, and `run_pipeline`
calls it once per genome. In a batch that is the same parse repeated.

Measured against the live database (101 records): 465 ms per call, against a
2.5 ms stat-only staleness stamp -- 186x cheaper to CHECK than to redo. Over
the 3,174 Basidiomycota genomes in BFD that is ~27 minutes of repeated YAML
parsing per rollout.

The cache must never serve a stale answer: the curated DB is edited by hand
between runs, and a detection run against a half-cached index would attribute
hits to the wrong family.
"""
import yaml

from MATPredict.detect import family_registry as fr
from MATPredict.detect.family_registry import FamilyKey, load_record_families


def _record(root, record_id, locus_name, phylum="Basidiomycota", order="Agaricales"):
    d = root / phylum / order / record_id
    d.mkdir(parents=True, exist_ok=True)
    (d / "metadata.yaml").write_text(
        yaml.safe_dump({"record_id": record_id, "mating_type": {"locus_name": locus_name}})
    )
    return d


def _counting_loads(monkeypatch) -> list[int]:
    calls: list[int] = []
    real = yaml.safe_load

    def counted(text):
        calls.append(1)
        return real(text)

    monkeypatch.setattr(fr.yaml, "safe_load", counted)
    return calls


def setup_function():
    fr.clear_record_families_cache()


def teardown_function():
    fr.clear_record_families_cache()


def test_repeated_calls_parse_the_database_once(tmp_path, monkeypatch):
    _record(tmp_path, "r1", "HD")
    _record(tmp_path, "r2", "PR")
    calls = _counting_loads(monkeypatch)
    first = load_record_families(tmp_path)
    for _ in range(4):
        assert load_record_families(tmp_path) == first
    assert len(calls) == 2  # two records, parsed once each


def test_the_answer_is_unchanged_by_caching(tmp_path):
    _record(tmp_path, "r1", "HD")
    _record(tmp_path, "r2", "PR")
    expected = {
        "r1": FamilyKey("Basidiomycota", "HD"),
        "r2": FamilyKey("Basidiomycota", "PR"),
    }
    assert load_record_families(tmp_path) == expected
    assert load_record_families(tmp_path) == expected


def test_an_edited_record_invalidates_the_cache(tmp_path, monkeypatch):
    """A curator changing a record's locus_name between runs must be seen."""
    _record(tmp_path, "r1", "HD")
    assert load_record_families(tmp_path)["r1"] == FamilyKey("Basidiomycota", "HD")
    _record(tmp_path, "r1", "PR")
    assert load_record_families(tmp_path)["r1"] == FamilyKey("Basidiomycota", "PR")


def test_a_new_record_invalidates_the_cache(tmp_path):
    _record(tmp_path, "r1", "HD")
    assert set(load_record_families(tmp_path)) == {"r1"}
    _record(tmp_path, "r2", "PR")
    assert set(load_record_families(tmp_path)) == {"r1", "r2"}


def test_a_removed_record_invalidates_the_cache(tmp_path):
    _record(tmp_path, "r1", "HD")
    d = _record(tmp_path, "r2", "PR")
    assert set(load_record_families(tmp_path)) == {"r1", "r2"}
    (d / "metadata.yaml").unlink()
    assert set(load_record_families(tmp_path)) == {"r1"}


def test_two_database_roots_are_cached_separately(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    _record(a, "r1", "HD")
    _record(b, "r2", "PR")
    assert set(load_record_families(a)) == {"r1"}
    assert set(load_record_families(b)) == {"r2"}
    assert set(load_record_families(a)) == {"r1"}
