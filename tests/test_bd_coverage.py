"""Tests for the Batrachochytrium depth, heterozygosity and parent-of-origin scripts."""
import gzip
import importlib.util
import pathlib

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


po = _load("bd_parent_of_origin")
bs = _load("bd_summarize")


def _paf(tmp_path, lines):
    p = tmp_path / "x.paf.gz"
    with gzip.open(p, "wt") as fo:
        fo.write("\n".join(lines) + "\n")
    return str(p)


def _line(tstart, cs, block=6000, tp="P", tname="T1"):
    return "\t".join(["Q", "9000", "0", "6000", "+", tname, "9000", str(tstart), str(tstart + 6000), "5990", str(block),
                      "60", f"tp:A:{tp}", f"cs:Z:{cs}"])


def test_divergent_sites_positions_and_deletions(tmp_path):
    # 10 matches, substitution c->t at 1-based 11, 5 matches, 2-base deletion (target only), 3 matches, sub a->g at 1-based 22
    paf = _paf(tmp_path, [_line(0, ":10*ct:5-aa:3*ag")])
    assert po.divergent_sites(paf, 5000) == {("T1", 11): "T", ("T1", 22): "G"}


def test_divergent_sites_filters(tmp_path):
    short = _line(0, ":10*ct", block=100)
    secondary = _line(0, ":10*ct", tp="S")
    assert po.divergent_sites(_paf(tmp_path, [short, secondary]), 5000) == {}


def test_divergent_site_covered_twice_is_dropped(tmp_path):
    paf = _paf(tmp_path, [_line(0, ":10*ct"), _line(0, ":10*ct")])
    assert po.divergent_sites(paf, 5000) == {}


def test_site_calls_absent_is_not_reference_allele():
    sites = {("T1", 11): "T", ("T1", 12): "G", ("T1", 13): "A"}
    snps = {("T1", 11): ("T", 5, 15), ("T1", 12): ("C", 4, 16)}  # 12 has another alt base
    depth = {("T1", 0): 30.0}
    out = {p: f for _, p, f in po.site_calls(sites, snps, depth, 10)}
    assert out == {11: 0.75, 13: None}


def test_site_calls_low_depth_window_skipped():
    assert list(po.site_calls({("T1", 11): "T"}, {}, {("T1", 0): 4.0}, 10)) == []


def test_window_table_classes():
    calls = [("T1", p, None) for p in range(1, 21)]  # 20 sites, none called -> CLFT
    calls += [("T1", 50000 + p, 1.0) for p in range(1, 21)]  # all hom JEL
    calls += [("T1", 100000 + p, 0.5) for p in range(1, 21)]  # all het
    calls += [("T1", 150000 + p, 1.0) for p in range(1, 6)]  # too few sites: dropped
    rows = po.window_table(calls, 50000, 20)
    assert [r[6] for r in rows] == ["CLFT", "JELHOM", "HET"]
    assert rows[1][4] == 1.0 and rows[2][5] == 1.0
    assert po.switches(rows) == 2


def test_switches_not_counted_across_sequences():
    rows = [("A", 0, 20, 20, 1.0, 0.0, "JELHOM"), ("B", 0, 20, 0, 0.0, 0.0, "CLFT")]
    assert po.switches(rows) == 0


def _write_run(base, ref, run, depths, snps):
    d = base / f"per_run_{ref}"
    d.mkdir(parents=True, exist_ok=True)
    with gzip.open(d / f"{run}.regions.bed.gz", "wt") as fo:
        for i, v in enumerate(depths):
            fo.write(f"S1\t{i * 10000}\t{(i + 1) * 10000}\t{v}\n")
    with gzip.open(d / f"{run}.snps.tsv.gz", "wt") as fo:
        for pos, ref_, alt, dp, ad in snps:
            fo.write(f"S1\t{pos}\t{ref_}\t{alt}\t50\t{dp}\t{ad}\n")


def test_summarise_run_het_and_depth(tmp_path):
    depths = [30.0] * 20 + [60.0] * 20  # second half doubled
    snps = [(100, "A", "G", 30, "15,15"), (200, "A", "G", 30, "1,29"), (300, "A", "C,G", 30, "10,10,10"),
            (400, "A", "G", 5, "2,3"), (500, "A", "G", 30, "30,0")]
    _write_run(tmp_path, "R", "RUN1", depths, snps)
    row, seqs, wrows = bs.summarise_run("RUN1", str(tmp_path), "R")
    # rows kept: 15,15 (het), 1,29 (hom alt), 30,0 (alt fraction 0); multiallelic and DP<10 rows are dropped
    assert row["snps"] == 3
    assert row["het_snps"] == 1
    assert row["hom_alt_fraction"] == round(1 / 3, 4)
    assert row["median_depth"] == 45.0
    assert len(wrows) == 40
    assert seqs[0][3] == 1.0  # one sequence: its median is the genome median
