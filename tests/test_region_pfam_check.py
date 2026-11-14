"""Tests for scripts/region_pfam_check.py (translation and ORF helpers) and scripts/matpi_tblastn.py (merge, distance)."""
import importlib.util
import pathlib

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


rp = _load("region_pfam_check")
mt = _load("matpi_tblastn")


def test_translate_and_revcomp():
    assert rp.translate("ATGGCCTAA") == "MA*"
    assert rp.revcomp("AACG") == "CGTT"


def test_orfs_forward_and_reverse_coordinates():
    prot = "M" + "A" * 49  # 50 aa
    dna = "ATG" + "GCC" * 49
    seq = "TT" + dna + "TAA" + "CC"
    fw = [o for o in rp.orfs(seq, 40) if o[0] > 0]
    assert any(len(o[3]) >= 50 and o[3].startswith("MAAA") for o in fw)
    rc = rp.revcomp(seq)
    rv = [o for o in rp.orfs(rc, 40) if o[0] < 0]
    hit = [o for o in rv if len(o[3]) >= 50 and o[3].startswith("MAAA")][0]
    # the ORF spans the same 150 bases of the reverse-complemented input: start0 .. end0 inside the sequence
    assert 0 <= hit[1] < hit[2] <= len(rc)
    assert hit[2] - hit[1] >= 150


def test_orfs_min_length_filter():
    assert rp.orfs("ATGGCCTAA", 40) == []


def test_merge_hits_and_distance():
    hits = [("c1", 100, 300, 50.0, 1e-10, 80.0), ("c1", 1000, 1200, 60.0, 1e-20, 90.0), ("c1", 9000, 9100, 40.0, 1e-5, 30.0),
            ("c2", 500, 400, 45.0, 1e-8, 70.0)]
    regs = mt.merge_hits(hits, 2000)
    assert [(r["contig"], r["start"], r["end"]) for r in regs] == [("c1", 100, 1200), ("c2", 400, 500), ("c1", 9000, 9100)]
    assert regs[0]["ident"] == 60.0 and regs[0]["bits"] == 170.0
    assert mt.distance(100, 200, 500, 600) == 300
    assert mt.distance(100, 700, 500, 600) == 0
