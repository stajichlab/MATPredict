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
    """A genome whose every 8-character block is unique (digits stand in for bases), so a window's start can be found by search."""
    g = tmp_path / "g.fa"
    seq = "".join(f"{i:08d}" for i in range(length // 8))
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
    assert S._tblastn_twostage_window_bp() == 60_000
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

def _runner(calls, core_hits, genome_flank_hits, seq=None):
    """tblastn stand-in. Core queries hit genome-wide. Flank queries: genome-wide when searched against the genome, and against a windows
    database only the genome-coordinate hits (`genome_flank_hits`: [(query, contig, start, end)]) that lie wholly inside one window,
    reported in window coordinates, as a real tblastn would."""
    def run(cmd, **kw):
        calls.append(cmd)
        if cmd[0] == "makeblastdb":
            return _res("")
        if cmd[0] == "blastdbcmd":
            return _res("Database: x\n\t2 sequences; 402,000 total bases\n")
        q = Path(cmd[cmd.index("-query") + 1]).read_text()
        if "MAT1-1-1" in q:
            return _res(core_hits)
        if "-dbsize" not in cmd:
            return _res("".join(_line(qq, c, s, e) for qq, c, s, e in genome_flank_hits))
        windows = Path(cmd[cmd.index("-db") + 1]).parent / "windows.fasta"
        out, name = "", None
        for line in windows.read_text().splitlines():
            if line.startswith(">"):
                name = line[1:]
            else:
                start = seq.find(line[:16]) + 1
                end = start + len(line) - 1
                for qq, c, s, e in genome_flank_hits:
                    if c == "chr1" and start <= s and e <= end:
                        out += _line(qq, name, s - start + 1, e - start + 1)
        return _res(out)
    return run


def _hits(tmp_path, runner):
    ref = tmp_path / "ref.faa"
    ref.write_text(REF)
    return search_localize(_genome(tmp_path), [FAM], ref, RECS, runner=runner)


def _seq(length=400_000):
    return "".join(f"{i:08d}" for i in range(length // 8))


def test_off_makes_exactly_one_genome_wide_call(tmp_path, monkeypatch):
    monkeypatch.delenv("MATPREDICT_TBLASTN_TWOSTAGE", raising=False)
    calls = []
    _hits(tmp_path, _runner(calls, _line("rec|gene0|MAT1-1-1", "chr1", 200_000, 200_300), [], _seq()))
    assert [c[0] for c in calls].count("tblastn") == 1 and not any("-dbsize" in c for c in calls)


def test_two_stage_searches_flanks_in_windows_and_maps_coordinates_back(tmp_path, monkeypatch):
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE", "yes")
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE_WINDOW_BP", "50000")
    calls = []
    flank = [("rec|gene1|SLA2", "chr1", 180_000, 180_300)]          # 20 kb from the core hit: inside the first window
    hits = _hits(tmp_path, _runner(calls, _line("rec|gene0|MAT1-1-1", "chr1", 200_000, 200_300), flank, _seq()))
    by = {h.gene_name: h for h in hits}
    assert (by["MAT1-1-1"].contig, by["MAT1-1-1"].start, by["MAT1-1-1"].end) == ("chr1", 200_000, 200_300)
    assert (by["SLA2"].contig, by["SLA2"].start, by["SLA2"].end) == ("chr1", 180_000, 180_300)
    tb = [c for c in calls if c[0] == "tblastn"]
    assert "-dbsize" not in tb[0] and all(c[c.index("-dbsize") + 1] == "402000" for c in tb[1:])


def test_two_stage_grows_the_windows_along_a_cluster_to_a_fixed_point(tmp_path, monkeypatch):
    """A flank gene 90 kb from the core hit is outside the 50 kb window, but a nearer flank gene (40 kb) widens the window to reach it,
    as the pipeline's gap rule would join them; a gene 200 kb away with nothing between is not found (it is not in the cluster)."""
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE", "yes")
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE_WINDOW_BP", "50000")
    calls = []
    flank = [("rec|gene1|SLA2", "chr1", 160_000, 160_300), ("rec|gene2|APN2", "chr1", 110_000, 110_300), ("rec|gene2|APN2", "chr1", 399_000, 399_300)]
    hits = _hits(tmp_path, _runner(calls, _line("rec|gene0|MAT1-1-1", "chr1", 200_000, 200_300), flank, _seq()))
    found = sorted((h.gene_name, h.start) for h in hits)
    assert found == [("APN2", 110_000), ("MAT1-1-1", 200_000), ("SLA2", 160_000)]
    assert len([c for c in calls if c[0] == "tblastn"]) >= 3                # core + at least two window rounds (growth, then no change)


def test_two_stage_without_a_core_hit_searches_the_flanks_genome_wide(tmp_path, monkeypatch):
    monkeypatch.setenv("MATPREDICT_TBLASTN_TWOSTAGE", "yes")
    calls = []
    hits = _hits(tmp_path, _runner(calls, "", [("rec|gene1|SLA2", "chr2", 100, 400)], _seq()))
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
    search_localize(_genome(tmp_path), [fam], ref, RECS, runner=_runner(calls, _line("rec|gene0|MAT1-1-1", "chr1", 5, 300), [], _seq()))
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
        _hits(tmp_path, _runner(calls, "", [], _seq()))
    except Exception:
        pass
    assert not any(c[0] == "blastdbcmd" for c in calls)


def test_the_options_page_names_the_setting():
    assert "MATPREDICT_TBLASTN_TWOSTAGE" in (Path(__file__).resolve().parents[2] / "docs" / "detect-options.md").read_text()
