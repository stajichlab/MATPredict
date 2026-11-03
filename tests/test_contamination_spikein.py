"""The contamination spike-in: mixes are random but reproducible, any prefix keeps the fraction, and the detection
logic counts a contaminant only when GmG names it with enough evidence."""
import gzip
import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("spikein", REPO / "datasets/reads/contamination_spikein/spikein.py")
sp = importlib.util.module_from_spec(spec)
sys.modules["spikein"] = sp
spec.loader.exec_module(sp)


def _fastq(path: Path, tag: str, n: int) -> tuple[str, str]:
    """n pairs whose read names carry `tag`, so a mixed file shows where each pair came from."""
    r1, r2 = path / f"{tag}_R1.fastq.gz", path / f"{tag}_R2.fastq.gz"
    for p, mate in ((r1, 1), (r2, 2)):
        with gzip.open(p, "wt") as fh:
            for i in range(n):
                fh.write(f"@{tag}{i}/{mate}\nACGTACGTAC\n+\nIIIIIIIIII\n")
    return str(r1), str(r2)


def _names(path: str) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return [line.strip() for i, line in enumerate(fh) if i % 4 == 0]


def test_the_mix_has_the_asked_fraction_and_keeps_each_sources_file_order(tmp_path):
    bg, con = _fastq(tmp_path, "bg", 20000), _fastq(tmp_path, "con", 20000)
    o1, o2 = str(tmp_path / "m_R1.fastq.gz"), str(tmp_path / "m_R2.fastq.gz")
    meta = sp.mix(bg, con, 0.1, 10000, seed=7, out_r1=o1, out_r2=o2, checkpoints=(1000, 10000))
    names = _names(o1)
    assert len(names) == 10000
    from_con = [n for n in names if n.startswith("@con")]
    assert meta["contaminant_pairs"] == len(from_con)
    assert abs(len(from_con) / 10000 - 0.1) < 0.02
    # each source's pairs appear in file order, with no gaps
    assert [n for n in names if n.startswith("@con")] == [f"@con{i}/1" for i in range(len(from_con))]
    assert [n for n in names if n.startswith("@bg")] == [f"@bg{i}/1" for i in range(10000 - len(from_con))]
    # mates stay together
    assert [n.replace("/2", "/1") for n in _names(o2)] == names


def test_any_prefix_is_a_mix_at_about_the_fraction_and_the_checkpoints_are_exact(tmp_path):
    bg, con = _fastq(tmp_path, "bg", 30000), _fastq(tmp_path, "con", 30000)
    o1, o2 = str(tmp_path / "m_R1.fastq.gz"), str(tmp_path / "m_R2.fastq.gz")
    meta = sp.mix(bg, con, 0.05, 20000, seed=3, out_r1=o1, out_r2=o2, checkpoints=(2000, 5000, 20000))
    names = _names(o1)
    for n, k in meta["checkpoints"].items():
        assert sum(1 for x in names[:n] if x.startswith("@con")) == k
        assert abs(k / n - 0.05) < 0.02
    assert set(meta["checkpoints"]) == {2000, 5000, 20000}


def test_the_same_seed_gives_the_same_mix_and_a_new_seed_a_different_one(tmp_path):
    bg, con = _fastq(tmp_path, "bg", 5000), _fastq(tmp_path, "con", 5000)

    def run(seed, tag):
        o1, o2 = str(tmp_path / f"{tag}_R1.fastq.gz"), str(tmp_path / f"{tag}_R2.fastq.gz")
        sp.mix(bg, con, 0.2, 3000, seed=seed, out_r1=o1, out_r2=o2)
        return _names(o1)
    assert run(1, "a") == run(1, "b")
    assert run(1, "a") != run(2, "c")


def test_a_zero_fraction_control_is_the_background_only(tmp_path):
    bg, con = _fastq(tmp_path, "bg", 1000), _fastq(tmp_path, "con", 1000)
    o1, o2 = str(tmp_path / "m_R1.fastq.gz"), str(tmp_path / "m_R2.fastq.gz")
    meta = sp.mix(bg, con, 0.0, 500, seed=1, out_r1=o1, out_r2=o2)
    assert meta["contaminant_pairs"] == 0
    assert all(n.startswith("@bg") for n in _names(o1))


def test_a_source_that_runs_out_stops_with_a_clear_message(tmp_path):
    bg, con = _fastq(tmp_path, "bg", 1000), _fastq(tmp_path, "con", 10)
    with pytest.raises(SystemExit, match="contaminant reads ran out"):
        sp.mix(bg, con, 0.5, 500, seed=1, out_r1=str(tmp_path / "x1.gz"), out_r2=str(tmp_path / "x2.gz"))


def test_conditions_cover_every_pair_and_fraction_plus_a_control_per_background():
    sources = [
        {"id": "pos", "role": "background", "species": "Coccidioides posadasii", "group": "Cocc"},
        {"id": "imm", "role": "background", "species": "Coccidioides immitis", "group": "Cocc"},
        {"id": "hist", "role": "contaminant", "species": "Histoplasma capsulatum", "group": "Hist"},
        {"id": "imm_c", "role": "contaminant", "species": "Coccidioides immitis", "group": "Cocc"},
    ]
    conds = sp.conditions(sources, fractions=(0.01, 0.1))
    assert len(conds) == 2 * (1 + 2 * 2)
    assert sum(1 for c in conds if c["kind"] == "control") == 2
    assert {c["kind"] for c in conds if c["con"] == "imm_c"} == {"same_genus"}
    assert {c["kind"] for c in conds if c["con"] == "hist"} == {"other_genus"}
    assert len({c["id"] for c in conds}) == len(conds) and len({c["seed"] for c in conds}) == len(conds)


def _gmg(rid, primary, second=None, hashes=0, ref_cov=0.0, suspected=False, conc=0.5):
    return {"id": rid, "status": "ok", "best_rank": "genus", "qc_flags": ["mixture_suspected"] if suspected else [],
            "guess": {"calls": {"genus": {"taxon": primary}, "species": {"taxon": None, "vote": primary + " sp."}}},
            "qc": {"coverage_est": 9.0, "n_hashes": 2500, "mixture": {"primary": primary, "explained_fraction": 0.9, "suspected": suspected,
                   "second": {"taxon": second, "hashes": hashes, "explained_fraction": hashes / 2500, "ref_coverage": ref_cov, "concentration": conc} if second else None}}}


def test_gmg_row_reads_the_nested_fields_and_tolerates_a_missing_second_genus():
    r = sp.gmg_row(_gmg("a", "Coccidioides", "Histoplasma", 300, 0.2, True))
    assert (r["genus"], r["second"], r["second_hashes"], r["second_ref_coverage"], r["mixture_suspected"]) == ("Coccidioides", "Histoplasma", 300, 0.2, True)
    r = sp.gmg_row(_gmg("b", "Coccidioides"))
    assert r["second"] is None and r["second_hashes"] is None and r["mixture_suspected"] is False


def test_sees_needs_the_contaminant_genus_with_enough_hashes_and_concentration():
    row = sp.gmg_row(_gmg("a", "Coccidioides", "Histoplasma", 60, 0.05, conc=0.9))
    assert sp.sees(row, "Histoplasma", 50, 0.8)
    assert not sp.sees(row, "Histoplasma", 100, 0.8)       # too few hashes
    assert not sp.sees(row, "Histoplasma", 50, 0.95)       # hashes not concentrated enough in one genome
    assert not sp.sees(row, "Lecanicillium", 5, 0.5)       # a different genus is not the contaminant
    noise = sp.gmg_row(_gmg("n", "Coccidioides", "Aspergillus", 12, 0.0008, conc=0.17))
    assert not sp.sees(noise, "Aspergillus", 5, 0.5)       # the scattered hashes of a clean set
    pilot = sp.gmg_row(_gmg("p", "Coccidioides", "Histoplasma", 6, 0.0016, conc=0.83))
    assert sp.sees(pilot, "Histoplasma", 5, 0.7) and not sp.sees(pilot, "Histoplasma", 20, 0.7)   # visible at 5%, 100k pairs, below GmG's own flag


def test_false_positives_count_any_second_genus_over_the_thresholds():
    neg = [sp.gmg_row(_gmg("c1", "Coccidioides", "Aspergillus", 12, 0.0008, conc=0.17)), sp.gmg_row(_gmg("c2", "Coccidioides", "Aspergillus", 60, 0.0015, conc=0.08)),
           sp.gmg_row(_gmg("c3", "Coccidioides"))]
    fp = sp.false_positives(neg, [(5, 0.5), (50, 0.0)])
    assert fp["n"] == 3 and fp["gmg_flag"] == 0
    assert fp["h5_q0.5"] == 0 and fp["h50_q0"] == 1


def test_a_contaminant_can_be_limited_to_some_backgrounds():
    sources = [
        {"id": "pos", "role": "background", "species": "Coccidioides posadasii", "group": "Cocc"},
        {"id": "imm", "role": "background", "species": "Coccidioides immitis", "group": "Cocc"},
        {"id": "imm_c", "role": "contaminant", "species": "Coccidioides immitis", "group": "Cocc", "only_with": "pos"},
    ]
    conds = sp.conditions(sources, fractions=(0.01,))
    assert [(c["bg"], c["con"]) for c in conds] == [("pos", "none"), ("pos", "imm_c"), ("imm", "none")]


def test_mix_one_builds_a_condition_from_the_sheet_and_is_idempotent(tmp_path):
    import csv as _csv
    bg, con = _fastq(tmp_path, "bg", 3000), _fastq(tmp_path, "con", 3000)
    src = tmp_path / "sources.tsv"
    src.write_text("id\tr1\tr2\tspecies\trole\tgroup\n" + f"b\t{bg[0]}\t{bg[1]}\tCoccidioides posadasii\tbackground\tC\n" + f"h\t{con[0]}\t{con[1]}\tHistoplasma capsulatum\tcontaminant\tH\n")
    out = tmp_path / "out"
    sp.main(["sheet", "--sources", str(src), "--outdir", str(out), "--fractions", "0.1"])
    rows = list(_csv.DictReader(open(out / "conditions.tsv"), delimiter="\t"))
    assert [r["id"] for r in rows] == ["b__none__f0", "b__h__f0.1"]
    sheet = (out / "samples.tsv").read_text().splitlines()
    assert sheet[0] == "id\tpath\texpected" and sheet[1].endswith("Coccidioides posadasii") and len(sheet) == 5   # two conditions, two files each
    sp.main(["mix-one", "--sources", str(src), "--outdir", str(out), "--index", "2", "--pairs", "2000"])
    meta = json.loads((out / "mixed" / "b__h__f0.1_R1.fastq.gz.mix.json").read_text()) if (out / "mixed" / "b__h__f0.1_R1.fastq.gz.mix.json").exists() else json.loads((out / "mixed" / "b__h__f0.1.mix.json").read_text())
    assert meta["id"] == "b__h__f0.1" and meta["n_pairs"] == 2000 and 100 < meta["contaminant_pairs"] < 320
    names = _names(str(out / "mixed" / "b__h__f0.1_R1.fastq.gz"))
    assert len(names) == 2000 and sum(1 for n in names if n.startswith("@con")) == meta["contaminant_pairs"]
    before = (out / "mixed" / "b__h__f0.1_R1.fastq.gz").stat().st_mtime_ns
    sp.main(["mix-one", "--sources", str(src), "--outdir", str(out), "--index", "2", "--pairs", "2000"])       # skipped
    assert (out / "mixed" / "b__h__f0.1_R1.fastq.gz").stat().st_mtime_ns == before
    sp.main(["mix-one", "--sources", str(src), "--outdir", str(out), "--index", "1", "--pairs", "1000"])       # a control
    assert json.loads((out / "mixed" / "b__none__f0.mix.json").read_text())["contaminant_pairs"] == 0
