"""datasets/reads/sra_survey/pilot/clone_solid.py: solid-containment ANI ignores sequencing errors and coverage gaps, average-linkage
clone labels, and the per-species tables (no sourmash needed: sketches are {hash: abundance} dicts)."""
import csv, importlib.util, random, sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
REPO = Path(__file__).resolve().parent.parent
PATH = REPO / "datasets" / "reads" / "sra_survey" / "pilot" / "clone_solid.py"
spec = importlib.util.spec_from_file_location("clone_solid", PATH); cs = importlib.util.module_from_spec(spec); spec.loader.exec_module(cs)


def genome(seed, n=4000):
    rnd = random.Random(seed)
    return [rnd.randrange(1, 2 ** 50) for _ in range(n)]


def sample(hashes, seed, cov=20, errors=400, gaps=0.0005):
    """A sketch of one sequenced genome: abundances around cov, `errors` error hashes seen once or twice, a fraction `gaps` not seen at all (at 20x almost none)."""
    rnd = random.Random(seed)
    d = {h: max(1, int(rnd.gauss(cov, cov ** 0.5))) for h in hashes if rnd.random() >= gaps}
    for _ in range(errors): d[rnd.randrange(2 ** 50, 2 ** 51)] = rnd.choice((1, 2))
    return d


def strain(base, seed, frac):
    """Another strain: a fraction frac of the base genome's hashes replaced."""
    rnd = random.Random(seed)
    return [h if rnd.random() >= frac else rnd.randrange(2 ** 51, 2 ** 52) for h in base]


def test_same_genome_sequenced_twice_is_near_one_and_a_different_strain_is_not():
    g = genome(1); other = strain(g, 2, 0.05)
    M = cs.ani_matrix([sample(g, 10), sample(g, 11), sample(other, 12)])
    assert M[0, 1] > 0.9995                       # error hashes do not count as differences
    assert 0.995 < M[0, 2] < 0.999                # 5% of hashes differ: 0.95^(1/21) ~ 0.9976
    assert np.allclose(M, M.T) and np.all(np.diag(M) == 1)


def test_solid_drops_low_abundance_hashes():
    med, sol = cs.solid({1: 20, 2: 22, 3: 1, 4: 2, 5: 21})
    assert med == 20 and sol == {1, 2, 5}


def test_downsample_keeps_hashes_below_the_scaled_limit():
    lim = cs.MAXH // 1000
    assert cs.downsample({lim - 1: 5, lim: 5, 3: 1}) == {lim - 1: 5, 3: 1}


def test_clusters_label_the_largest_group_c1_and_cut_at_each_threshold():
    g = genome(3); far = strain(g, 4, 0.30); near = strain(g, 5, 0.03)
    sk = [sample(far, 20), sample(g, 21), sample(g, 22), sample(g, 23), sample(near, 24)]
    M = cs.ani_matrix(sk); labs = cs.clusters(M, (0.9995, 0.99))
    assert labs[0.9995][1] == labs[0.9995][2] == labs[0.9995][3] == "c1"
    assert len(set(labs[0.9995])) == 3                       # the clone of three, the near strain alone, the far strain alone
    assert labs[0.99][4] == "c1" and labs[0.99][0] != "c1"   # at 0.99 the near strain (0.97^(1/21) ~ 0.9985) joins
    mins = cs.group_min(M, labs[0.9995])
    assert mins[0] is None and mins[1] == mins[2] > 0.9995


def test_one_run_is_its_own_clone():
    assert cs.clusters(np.ones((1, 1)), (0.999,)) == {0.999: ["c1"]}


def test_species_groups_reads_every_committed_wave_once(tmp_path):
    for w, rows in (("wave1", [("SRR1", "sp_a"), ("SRR2", "sp_b")]), ("wave3_asco", [("SRR3", "sp_a"), ("SRR1", "sp_a")]), ("waveX_old", [("SRR9", "sp_a")])):
        (tmp_path / w).mkdir()
        with open(tmp_path / w / "results_all_final.tsv", "w") as fh:
            fh.write("run\tgroup\n" + "".join(f"{r}\t{g}\n" for r, g in rows))
    assert dict(cs.species_groups(str(tmp_path))) == {"sp_a": ["SRR1", "SRR3"], "sp_b": ["SRR2"]}


def test_write_group_gives_a_row_per_run_and_a_summary(tmp_path):
    g = genome(6)
    s = cs.write_group("sp_a", ["R1", "R2", "R3"], [sample(g, 1), sample(g, 2), sample(strain(g, 7, 0.1), 3, cov=5, gaps=0.01)], str(tmp_path))
    rows = list(csv.DictReader(open(tmp_path / "sp_a.clones.tsv"), delimiter="\t"))
    assert [r["run"] for r in rows] == ["R1", "R2", "R3"]
    assert rows[0]["nearest_run"] == "R2" and rows[0]["clone_0.9995"] == rows[1]["clone_0.9995"] == "c1"
    assert rows[2]["min_ani_0.9995"] == "" and float(rows[2]["median_abund"]) < cs.LOW_COVERAGE
    assert s["runs"] == 3 and s["clones_0.9995"] == 2 and s["low_coverage_runs"] == 1
    assert s["same_strain_pairs"] == 0


def test_same_strain_pairs_use_normalised_names_and_skip_placeholders():
    assert cs.strain_key("CBS 144.89") == cs.strain_key("cbs144-89") == "cbs14489"
    assert cs.strain_key("missing") == cs.strain_key("NA") == cs.strain_key("") == ""
    M = np.array([[1, .9999, .99], [.9999, 1, .98], [.99, .98, 1]])
    assert cs.same_strain(M, ["a", "b", "c"], {"a": "Af293", "b": "AF-293", "c": "missing"}) == [pytest.approx(.9999)]
