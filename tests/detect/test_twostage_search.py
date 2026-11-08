"""The opt-in two-stage tblastn search (`MATPREDICT_TBLASTN_TWOSTAGE`, #287 step 2)."""
from __future__ import annotations

from pathlib import Path

import pytest

from MATPredict.detect import search as S
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.search import search_localize

KEY = FamilyKey("Ascomycota", "MAT")
FAM = Family(key=KEY, vocabulary_type="enum", idiomorph_values=["MAT1-1", "MAT1-2"], idiomorph_pattern=None,
             genes=[{"name": "MAT1-1-1", "role": "core_MAT"}, {"name": "SLA2", "role": "flanking_conserved"},
                    {"name": "APN2", "role": "flanking_conserved"}], taxonomic_scope=[4890])
RECS = {"rec": KEY}
REF = ">rec|gene0|MAT1-1-1\nMKKKK\n>rec|gene1|SLA2\nMKKKKKK\n>rec|gene2|APN2\nMKKK\n"


def _res(stdout: str):
    class R:
        returncode = 0
        stderr = ""
    R.stdout = stdout
    return R()


def _genome(tmp_path, length=400_000):
    g = tmp_path / "g.fa"
    seq = ("ACGT" * (length // 4))
    g.write_text(f">chr1\n{seq}\n>chr2\n{seq[:2000]}\n")
    return g


def _line(q, contig, s, e, bits="80.0"):
    return f"{q}\t{contig}\t90.0\t100\t200\t50\t1e-30\t{bits}\t{s}\t{e}\t1\n"


# ---- the setting -------------------------------------------------------------------------------------------------------------------

def test_off_by_default_and_validated(monkeypatch):
    monkeypatch.delenv("MATPREDICT_TBLASTN_TWOSTAGE", raising=False)
    assert S._tblastn_twostage() is False
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE", "yes")
    assert S._tblastn_twostage() is True
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE", "maybe")
    with pytest.raises(ValueError):
        S._tblastn_twostage()


def test_window_size_default_and_validation(monkeypatch):
    monkeypatch.delenv("MATPREDICT_TBLASTN_TWOSTAGE_WINDOW_BP", raising=False)
    assert S._tblastn_twostage_window_bp() == 100_000
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE_WINDOW_BP", "5000")
    assert S._tblastn_twostage_window_bp() == 5000
    for bad in ("0", "-5", "ten"):
        monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE_WINDOW_BP", bad)
        with pytest.raises(ValueError):
            S._tblastn_twostage_window_bp()


# ---- helpers -----------------------------------------------------------------------------------------------------------------------

def test_split_by_role_keeps_the_file_order(tmp_path):
    ref = tmp_path / "ref.faa"
    ref.write_text(REF)
    core, rest, order = S._split_reference_by_role(ref, tmp_path, lambda q: q.endswith("MAT1-1-1"))
    assert order == ["rec|gene0|MAT1-1-1", "rec|gene1|SLA2", "rec|gene2|APN2"]
    assert core.read_text() == ">rec|gene0|MAT1-1-1\nMKKKK\n"
    assert rest.read_text() == ">rec|gene1|SLA2\nMKKKKKK\n>rec|gene2|APN2\nMKKK\n"


def test_split_by_role_returns_none_for_an_empty_side(tmp_path):
    ref = tmp_path / "ref.faa"
    ref.write_text(REF)
    core, rest, _ = S._split_reference_by_role(ref, tmp_path, lambda q: True)
    assert core is not None and rest is None


def test_windows_are_padded_clipped_and_merged():
    spans = [("c", 1000, 2000), ("c", 5000, 6000), ("c", 50_000, 51_000), ("d", 10, 20)]
    assert S._merge_windows(spans, 3000) == [("c", 1, 9000), ("c", 47_000, 54_000), ("d", 1, 3020)]


def test_overlapping_windows_on_one_contig_merge():
    assert S._merge_windows([("c", 100, 200), ("c", 250, 300)], 100) == [("c", 1, 400)]


def test_window_lines_go_back_to_genome_coordinates():
    out = "q\twin0\t90.0\t100\t200\t50\t1e-30\t80.0\t11\t310\t1\nq\twin1\t90.0\t100\t200\t50\t1e-30\t80.0\t400\t5\t-1\n"
    back = S._remap_window_lines(out, {"win0": ("chr1", 1001), "win1": ("chr2", 1)}).splitlines()
    assert back[0].split("\t")[1:] == ["chr1", "90.0", "100", "200", "50", "1e-30", "80.0", "1011", "1310", "1"][0:1] + \
        "90.0 100 200 50 1e-30 80.0 1011 1310 1".split()
    assert back[1].split("\t")[1] == "chr2" and back[1].split("\t")[8:10] == ["400", "5"]


# ---- end to end with a fake tblastn --------------------------------------------------------------------------------------------------

def _runner(calls, core_hits, flank_hits_in_windows, genome_flank_hits=""):
    """tblastn stand-in: core queries hit genome-wide; flank queries hit genome-wide, or in the window database when `-dbsize` is set."""
    def run(cmd, **kw):
        calls.append(cmd)
        if cmd[0] == "makeblastdb":
            return _res("")
        if cmd[0] == "blastdbcmd":
            return _res("Database: x\n\t2 sequences; 402,000 total bases\n")
        q = Path(cmd[cmd.index("-query") + 1]).read_text()
        is_core = "MAT1-1-1" in q
        if is_core:
            return _res(core_hits)
        return _res(flank_hits_in_windows if "-dbsize" in cmd else genome_flank_hits)
    return run


def _hits(tmp_path, runner):
    ref = tmp_path / "ref.faa"
    ref.write_text(REF)
    return search_localize(_genome(tmp_path), [FAM], ref, RECS, runner=runner)


def test_off_makes_exactly_one_genome_wide_call(tmp_path, monkeypatch):
    monkeypatch.delenv("MATPREDICT_TBLASTN_TWOSTAGE", raising=False)
    calls = []
    _hits(tmp_path, _runner(calls, _line("rec|gene0|MAT1-1-1", "chr1", 200_000, 200_300), ""))
    assert [c[0] for c in calls].count("tblastn") == 1 and not any("-dbsize" in c for c in calls)


def test_two_stage_searches_flanks_in_windows_and_maps_coordinates_back(tmp_path, monkeypatch):
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE", "yes")
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE_WINDOW_BP", "50000")
    calls = []
    # window 0 covers chr1 150000..250300 (core hit 200000..200300 padded by 50000); a flank hit at window offset 1000..1300 is chr1 150999..151299
    flank = _line("rec|gene1|SLA2", "win0", 1000, 1300)
    hits = _hits(tmp_path, _runner(calls, _line("rec|gene0|MAT1-1-1", "chr1", 200_000, 200_300), flank))
    by = {h.gene_name: h for h in hits}
    assert (by["MAT1-1-1"].contig, by["MAT1-1-1"].start, by["MAT1-1-1"].end) == ("chr1", 200_000, 200_300)
    assert (by["SLA2"].contig, by["SLA2"].start, by["SLA2"].end) == ("chr1", 150_000 + 999, 150_000 + 1299)
    tb = [c for c in calls if c[0] == "tblastn"]
    assert len(tb) == 2 and "-dbsize" in tb[1] and tb[1][tb[1].index("-dbsize") + 1] == "402000"
    assert "-dbsize" not in tb[0]


def test_two_stage_without_a_core_hit_searches_the_flanks_genome_wide(tmp_path, monkeypatch):
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE", "yes")
    calls = []
    hits = _hits(tmp_path, _runner(calls, "", "", genome_flank_hits=_line("rec|gene1|SLA2", "chr2", 100, 400)))
    assert [h.gene_name for h in hits] == ["SLA2"] and (hits[0].contig, hits[0].start) == ("chr2", 100)
    tb = [c for c in calls if c[0] == "tblastn"]
    assert len(tb) == 2 and not any("-dbsize" in c for c in tb)


def test_two_stage_falls_back_to_the_plain_search_when_a_side_is_empty(tmp_path, monkeypatch):
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE", "yes")
    fam = Family(key=KEY, vocabulary_type="enum", idiomorph_values=["MAT1-1"], idiomorph_pattern=None,
                 genes=[{"name": "MAT1-1-1", "role": "core_MAT"}], taxonomic_scope=[4890])
    ref = tmp_path / "ref.faa"
    ref.write_text(">rec|gene0|MAT1-1-1\nMKKKK\n")
    calls = []
    search_localize(_genome(tmp_path), [fam], ref, RECS, runner=_runner(calls, _line("rec|gene0|MAT1-1-1", "chr1", 5, 300), ""))
    assert [c[0] for c in calls].count("tblastn") == 1


def test_the_output_order_matches_a_single_call(tmp_path, monkeypatch):
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE", "yes")
    out = S._merge_tblastn_by_query_order(
        [_line("rec|gene0|MAT1-1-1", "chr1", 1, 9), _line("rec|gene2|APN2", "chr1", 5, 99) + _line("rec|gene1|SLA2", "chr1", 7, 77)],
        ["rec|gene0|MAT1-1-1", "rec|gene1|SLA2", "rec|gene2|APN2"])
    assert [l.split("\t")[0] for l in out.splitlines()] == ["rec|gene0|MAT1-1-1", "rec|gene1|SLA2", "rec|gene2|APN2"]


def test_two_stage_is_ignored_when_seg_filtering_is_on(tmp_path, monkeypatch):
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE", "yes")
    monkeypatch.setenv("MATPREDICT_TBLASTN_SEG", "yes")
    calls = []
    try:
        _hits(tmp_path, _runner(calls, "", ""))
    except Exception:
        pass
    assert not any(c[0] == "blastdbcmd" for c in calls)


def test_the_options_page_names_the_setting():
    assert "MATPREDICT_TBLASTN_TWOSTAGE" in (Path(__file__).resolve().parents[2] / "docs" / "detect-options.md").read_text()
