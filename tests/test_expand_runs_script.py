"""expand_runs.pick: round-robin across studies, nearest 3 Gb first, a per-study cap and an overall need."""
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("expand_runs", Path(__file__).resolve().parents[1] / "datasets/reads/availability/expand_runs.py")
er = importlib.util.module_from_spec(spec)
spec.loader.exec_module(er)


def _rows(study, n, gb=3.0):
    return [{"study_accession": study, "run_accession": f"{study}_{i}", "base_count": str(int((gb + i * 0.1) * 1e9))} for i in range(n)]


def test_round_robin_spreads_studies_and_respects_both_caps():
    rows = _rows("BIG", 40) + _rows("MID", 5) + _rows("ONE", 1)
    got = er.pick(rows, need=12, per_study=25)
    assert [r["study_accession"] for r in got[:3]] == ["BIG", "MID", "ONE"]       # one from each study before a second from any
    assert len(got) == 12
    got = er.pick(rows, need=100, per_study=10)
    assert sum(r["study_accession"] == "BIG" for r in got) == 10 and len(got) == 16    # BIG stops at the per-study cap
    assert er.pick(rows, need=0, per_study=10) == []


def test_nearest_three_gigabases_first_inside_a_study():
    rows = [{"study_accession": "S", "run_accession": "far", "base_count": "9000000000"}, {"study_accession": "S", "run_accession": "near", "base_count": "3100000000"}]
    assert er.pick(rows, need=1, per_study=5)[0]["run_accession"] == "near"


def test_species_file_limits_the_species_and_skips_comments(tmp_path):
    f = tmp_path / "sp.txt"
    f.write_text("# thin-evidence species\nMadurella mycetomatis\tsingle_locus\n\nFusarium sambucinum\n")
    only = er.read_species_file(str(f))
    assert only == {"Madurella mycetomatis", "Fusarium sambucinum"} and er.read_species_file(None) is None
    rows = [{"species": "Madurella mycetomatis", "action": "unvalidated"}, {"species": "Fusarium sambucinum", "action": "held"},
            {"species": "Botrytis cinerea", "action": "expand"}]
    assert list(er.select_untyped(rows, ["expand", "unvalidated"], only)) == ["Madurella mycetomatis"]     # held is not an action asked for; Botrytis is not listed
    assert sorted(er.select_untyped(rows, ["expand", "unvalidated"])) == ["Botrytis cinerea", "Madurella mycetomatis"]
