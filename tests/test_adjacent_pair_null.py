"""Tests for scripts/adjacent_pair_null.py (pair helpers and the p-value logic on a toy lineage)."""
import importlib.util
import pathlib

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
spec = importlib.util.spec_from_file_location("adjacent_pair_null", SCRIPTS / "adjacent_pair_null.py")
apn = importlib.util.module_from_spec(spec)
spec.loader.exec_module(apn)


def test_has_pair_window_and_contig():
    assert apn.has_pair([("c1", 10)], [("c1", 13)], 5, False)
    assert not apn.has_pair([("c1", 10)], [("c1", 17)], 5, False)
    assert not apn.has_pair([("c1", 10)], [("c2", 11)], 5, False)


def test_has_pair_same_gene_is_not_a_pair():
    assert not apn.has_pair([("c1", 10)], [("c1", 10)], 5, True)
    assert apn.has_pair([("c1", 10), ("c1", 12)], [("c1", 10), ("c1", 12)], 5, True)


def test_pair_strands():
    a = [("c1", 10, "+")]
    b = [("c1", 12, "+"), ("c1", 14, "-"), ("c1", 30, "+")]
    assert apn.pair_strands(a, b, 5) == ["same", "opposite"]


def _write(tmp_path, name, genes, cands):
    with open(tmp_path / f"{name}.genes_all.tsv", "w") as fo:
        fo.write("gene\tcontig\tstart\tend\tstrand\tindex\n")
        for i, g in enumerate(genes):
            fo.write(f"{g}\tc1\t{i * 10}\t{i * 10 + 5}\t+\t{i}\n")
    with open(tmp_path / f"{name}.candidates.tsv", "w") as fo:
        fo.write("gene\tcontig\tindex\tpfam\tpfam_evalue\tcontrol_class\tcontrol_id\tcontrol_evalue\n")
        for g, d in cands.items():
            fo.write(f"{g}\tc1\t0\t{d}\t1e-10\t\t\t\n")


def test_main_counts_observed_pair(tmp_path, monkeypatch, capsys):
    genes = [f"g{i}" for i in range(200)]
    _write(tmp_path, "G1", genes, {"g10": "HMG_box", "g12": "Homeodomain"})
    (tmp_path / "list.txt").write_text("G1\n")
    out = tmp_path / "out.tsv"
    monkeypatch.setattr("sys.argv", ["x", "--stage1", str(tmp_path), "--out", str(out), "--perms", "50", str(tmp_path / "list.txt") + ":L"])
    apn.main()
    rows = {l.split("\t")[1]: l.split("\t") for l in out.read_text().splitlines()[1:]}
    assert rows["HMG+Homeo"][3] == "1"  # observed genomes
    assert float(rows["HMG+Homeo"][4]) < 0.5  # a random pair in 200 genes is rare
    assert rows["2xHomeo"][3] == "0"
