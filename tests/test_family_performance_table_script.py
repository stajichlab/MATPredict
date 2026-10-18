"""scripts/family_performance_table.py: per-family counts from detection reports, and held-out scores."""
import importlib.util
import json
from pathlib import Path

import yaml

_ROOT = Path(__file__).resolve().parents[1]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


fpt = _load("family_performance_table")


def gene(name, start, end, strand, contig="c1", bits=100):
    return dict(gene=name, contig=contig, start=start, end=end, strand=strand, bitscore=bits)


def locus(family, idio, conf="high", cls="mat_locus", genes=(), contig="c1"):
    return dict(family=family, contig=contig, start=1, end=5000, idiomorph=idio, confidence=conf, locus_class=cls,
                gene_evidence=list(genes))


def report(root, rel, attempted, detected=(), **extra):
    d = root / rel
    d.mkdir(parents=True)
    (d / "detection_report.yaml").write_text(yaml.safe_dump(
        dict(families_attempted=list(attempted), detected=list(detected), **extra)))


RED = "Basidiomycota:redHD"
PR = "Basidiomycota:redPR"


def build(tmp_path):
    hd = lambda a, b, s1, s2, c="c1": [gene("HD1", *a, s1, c), gene("HD2", *b, s2, c)]
    # g1: divergent pair, 400 bp apart; g2: tandem, called twice (two loci); g3: undetermined partial; g4: not called
    report(tmp_path, "b1/5286_GCA_1", [RED, PR], [
        locus(RED, "B41", genes=hd((1000, 1900), (2301, 3200), "-", "+")), locus(PR, "A1")])
    report(tmp_path, "b1/5535_GCA_2", [RED, PR], [
        locus(RED, "B63", "medium", genes=hd((1000, 1900), (2000, 2900), "+", "+")),
        locus(RED, "B8", "low", genes=hd((9000, 9900), (9500, 9990), "-", "+", "c2"), contig="c2")])
    report(tmp_path, "b2/5537_GCA_3", [RED, PR], [
        locus(RED, "undetermined", "low", "partial_locus", genes=[gene("HD1", 1, 900, "-")])])
    report(tmp_path, "b2/5005_GCA_4", [RED, PR])
    report(tmp_path, "b2/4890_GCA_5", [])  # outside every scope: not_searched
    (tmp_path / "b2" / "9999_GCA_6").mkdir()  # a failed run: directory, no report
    return tmp_path


def test_counts_per_family(tmp_path):
    res = fpt.summarize_reports(build(tmp_path))
    assert res["genomes"] == 5 and res["no_report"] == ["9999_GCA_6"]
    s = res["stats"][RED].row()
    assert (s["in_scope"], s["called"], s["loci"]) == (4, 3, 4)
    assert (s["conf_high"], s["conf_medium"], s["conf_low"]) == (1, 1, 2)
    assert (s["best_high"], s["best_medium"], s["best_low"]) == (1, 1, 1)
    assert s["undetermined"] == 1 and s["partial"] == 1 and s["multi_locus"] == 1 and s["both_alleles"] == 1
    assert "B41=1" in s["idiomorph"] and "undetermined=1" in s["idiomorph"]
    p = res["stats"][PR].row()
    assert (p["in_scope"], p["called"]) == (4, 1)


def test_pair_arrangement(tmp_path):
    arr = fpt.summarize_reports(build(tmp_path))["stats"][RED].row()["arrangement"]
    assert "divergent=2" in arr and "tandem=1" in arr and "pair_incomplete=1" in arr
    # gaps: 400 (g1), 99 (g2 tandem, 2000-1900-1), 0 (the c2 pair overlaps)
    assert "gap_bp_min/median/max=0/99/400" in arr
    assert fpt.orientation(gene("HD1", 1, 10, "+"), gene("HD2", 20, 30, "-")) == ("convergent", 9)
    assert fpt.orientation(gene("HD1", 20, 30, "+"), gene("HD2", 1, 10, "+")) == ("tandem", 9)


def test_orientation_is_by_position_not_by_gene_name():
    # HD2 on the left, minus strand; HD1 on the right, plus strand: divergent whichever gene is which
    assert fpt.orientation(gene("HD1", 500, 900, "+"), gene("HD2", 1, 100, "-")) == ("divergent", 399)


def test_family_filter_and_repeats(tmp_path):
    root = build(tmp_path)
    report(root, "b3/5286_GCA_1", [RED], [locus(RED, "B1")])  # same genome name in another batch
    res = fpt.summarize_reports(root, {PR})
    assert res["repeats"] == 1 and set(res["stats"]) == {PR}


def test_cli_writes_tsv_and_md(tmp_path, capsys):
    root = build(tmp_path / "runs")
    out = tmp_path / "t.tsv"
    assert fpt.main([str(root), "--out", str(out), "--md", str(tmp_path / "t.md")]) == 0
    lines = out.read_text().splitlines()
    assert lines[0].split("\t")[:4] == ["family", "in_scope", "called", "loci"]
    assert len(lines) == 3 and "5 reports" in capsys.readouterr().out
    assert (tmp_path / "t.md").read_text().startswith("5 reports")


def _row(rid, radius, family, expected, labels, status="ok", enum=True, left=3, idio_left=3, loci=None):
    return dict(record_id=rid, radius=radius, phylum="Basidiomycota", subphylum="", family=family,
                system="heterothallic", expected=expected, withheld=1, status=status, enum_family=enum,
                pipeline_labels=labels, loci=[{}] * len(labels) if loci is None else loci,
                ref_family_left=left, ref_idiomorph_left=idio_left, found_at_truth=None)


ROWS = [
    _row("a", "genus", "Basidiomycota:redPR", "A1", ["A1"]),
    _row("b", "genus", "Basidiomycota:redPR", "A1", ["A2"]),
    _row("c", "genus", "Basidiomycota:redPR", "A1", ["undetermined"]),
    _row("d", "genus", "Basidiomycota:redPR", "A1", []),
    _row("e", "family", "Basidiomycota:redPR", "A1", [], left=0, idio_left=0),
    _row("f", "genus", "Basidiomycota:redPR", "A1", [], idio_left=0),
    _row("g", "order", "Basidiomycota:rustHD", "b1", ["undetermined"], enum=False),
    _row("h", "order", "Basidiomycota:rustHD", "b1", [], enum=False),
    _row("i", "genus", "Basidiomycota:redHD", "B1", [], status="detect_error"),
    _row("j", "genus", "Basidiomycota:redHD", "B1", [], status="no_reference_family", left=0),
]


def test_holdout_buckets():
    out = {(r["family"], r["radius"]): r for r in fpt.summarize_holdout([dict(r) for r in ROWS])}
    pr = out[("Basidiomycota:redPR", "genus")]
    assert (pr["rows"], pr["tested"], pr["correct"], pr["wrong"], pr["undetermined"], pr["no_locus"],
            pr["no_reference"]) == (5, 4, 1, 1, 1, 1, 1)
    assert out[("Basidiomycota:redPR", "family")]["no_reference"] == 1
    assert out[("Basidiomycota:rustHD", "order")]["correct"] == 1  # pattern family: locus found
    assert out[("Basidiomycota:rustHD", "order")]["no_locus"] == 1
    hd = out[("Basidiomycota:redHD", "genus")]
    assert hd["error"] == 1 and hd["no_reference"] == 1 and hd["tested"] == 0


def test_holdout_manifest_recomputes_references(tmp_path):
    rows = tmp_path / "rows"
    rows.mkdir()
    (rows / "a.json").write_text(json.dumps(_row("a", "genus", "Basidiomycota:redPR", "A1", [], left=5, idio_left=5)))
    man = tmp_path / "m.tsv"
    man.write_text("record_id\tgenus_same_family_left\tgenus_same_idiomorph_left\na\t4\t0\n")
    (r,) = fpt.load_holdout([rows], man)
    assert (r["ref_family_left"], r["ref_idiomorph_left"]) == (4, 0)
    assert fpt.summarize_holdout([r])[0]["no_reference"] == 1  # own idiomorph has no reference left


def test_verdicts_match_run_asco_holdout_summarize(tmp_path, capsys):
    """Parity with the runner's own summary, so the two cannot drift."""
    runner = _load("run_asco_holdout")
    paths = []
    for i, r in enumerate(ROWS):
        p = tmp_path / f"{i}.json"
        p.write_text(json.dumps(r))
        paths.append(str(p))
    runner.summarize(type("A", (), {"rows": paths, "manifest": None})())
    printed = capsys.readouterr().out.splitlines()
    for r in ROWS:
        line = next(l for l in printed if l.startswith(r["record_id"][:34].ljust(34) + " " + r["radius"]))
        assert fpt.holdout_verdict(dict(r)) in line.split()


def test_holdout_listing_is_limited_to_searched_families_and_shows_manifest_gaps(tmp_path):
    man = tmp_path / "m.tsv"
    man.write_text(
        "record_id\tphylum\tfamily\ttestable\tgenus_same_family_left\tgenus_same_idiomorph_left\n"
        "a\tBasidiomycota\tredPR\tTrue\t4\t2\n"
        "b\tBasidiomycota\tredPR\tTrue\t0\t0\n"
        "c\tBasidiomycota\tredPR\tFalse\t4\t2\n"
        "d\tBasidiomycota\twallMAT\tTrue\t0\t0\n")
    rows = [_row("a", "genus", "Basidiomycota:redPR", "A1", ["A1"], left=4, idio_left=2)]
    out = fpt.summarize_holdout(rows, fpt.manifest_gaps(man), {"Basidiomycota:redPR"})
    assert [(r["family"], r["radius"]) for r in out] == [("Basidiomycota:redPR", "genus")]
    r = out[0]
    assert (r["records"], r["no_assembly"], r["no_reference_not_run"], r["correct"]) == (3, 1, 1, 1)
    # a radius with no row at all still appears, with its gaps
    out = fpt.summarize_holdout([], fpt.manifest_gaps(man), {"Basidiomycota:wallMAT"})
    assert out[0]["rows"] == 0 and out[0]["no_reference_not_run"] == 1
