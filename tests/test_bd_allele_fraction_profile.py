"""Tests for scripts/bd_allele_fraction_profile.py."""
import gzip
import importlib.util
import os

_p = os.path.join(os.path.dirname(__file__), "..", "scripts", "bd_allele_fraction_profile.py")
_s = importlib.util.spec_from_file_location("bd_allele_fraction_profile", _p)
m = importlib.util.module_from_spec(_s)
_s.loader.exec_module(m)


def test_binom_mid_and_profile(tmp_path):
    assert abs(m.binom_mid(10, 0.0, 1.0) - 1.0) < 1e-9
    assert m.binom_mid(2) == 0.5  # k=1 only: f = 0.5
    sites = tmp_path / "s.tsv"
    sites.write_text("c1\t10\tA\nc1\t20\tA\nc1\t30\tA\nc1\t40\tA\n")
    cnt = tmp_path / "c.tsv.gz"
    with gzip.open(cnt, "wt") as fh:
        fh.write("contig\tpos\tn_clft\tn_jel\nc1\t10\t10\t0\nc1\t20\t5\t5\nc1\t30\t0\t12\nc1\t40\t3\t1\nc1\t99\t5\t5\n")
    ks = m.load_sites(str(sites))
    hist, n, mid, exp, win = m.profile(str(cnt), ks, 8, 100000)
    assert n == 3  # site 40 has depth 4 (< 8); site 99 is not in the list
    assert hist[0] == 1 and hist[5] == 1 and hist[9] == 1
    assert mid == 1
    assert len(win[("c1", 0)]) == 3
