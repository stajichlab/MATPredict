"""Polish selection with an identity tier above the cap.

Measured on the full Dothideomycete run (analysis/2026-10-05_dothideomycetes-full-run.md): the plain
rank (distinct genes, then identity, then hits) put a true locus hitting 2 genes at ~100% behind six
noise regions hitting 3 genes at 33-43%, so with a cap of 6 it was never polished (15 genomes, including
the IPO323 reference). The tier polishes every cluster at or above `strong_identity`, then fills the
rest of the cap in the plain rank. Replay over 2,722 genomes: 0 true loci lost, same polishing work.
"""
import types

from MATPredict.detect.pipeline import (
    DEFAULT_POLISH_STRONG_IDENTITY,
    rank_for_polish,
    run_pipeline,
    select_polish_clusters,
)

from tests.detect.test_cap_v3_strong_fragment import KEY, _cluster, _hit
from tests.detect.test_pipeline import _model, _tblastn, _write_order, _write_record
from tests.detect.test_polish_cluster_cap import ORDER

OTHER = ("P", "OTHER")


def _scenario():
    """Six noise clusters (3 genes at 40%) and one true locus (2 genes at 100%)."""
    noise = [_cluster(3, 40.0, f"n{i}") for i in range(6)]
    true = _cluster(2, 100.0, "true")
    return noise, true


def test_the_default_tier_is_fifty_percent():
    assert DEFAULT_POLISH_STRONG_IDENTITY == 50.0


def test_the_plain_cap_skips_a_strong_cluster_that_hits_few_genes():
    noise, true = _scenario()
    selected = select_polish_clusters({KEY: noise + [true]}, 6, {})
    assert (id(true), KEY) not in selected          # the failure the tier fixes
    assert len(selected) == 6


def test_the_tier_polishes_the_strong_cluster_and_still_fills_the_cap():
    noise, true = _scenario()
    selected = select_polish_clusters({KEY: noise + [true]}, 6, {}, strong_identity=50.0)
    assert (id(true), KEY) in selected
    assert len(selected) == 6                       # same work as the plain cap
    # the five slots left are the plain rank of the rest
    rest = rank_for_polish(noise, KEY, {})[:5]
    assert {(id(c), KEY) for c in rest} <= selected


def test_more_strong_clusters_than_the_cap_are_all_polished_and_nothing_else():
    strong = [_cluster(1, 90.0, f"s{i}") for i in range(8)]
    weak = [_cluster(3, 35.0, f"w{i}") for i in range(4)]
    selected = select_polish_clusters({KEY: strong + weak}, 6, {}, strong_identity=50.0)
    assert selected == {(id(c), KEY) for c in strong}
    assert len(selected) == 8                       # the one case that exceeds the cap


def test_with_no_strong_cluster_the_selection_is_the_plain_cap():
    weak = [_cluster(n, 35.0 + n, f"w{i}") for i, n in enumerate([1, 2, 3, 2, 1, 3, 2, 1])]
    assert (select_polish_clusters({KEY: weak}, 6, {}, strong_identity=50.0)
            == select_polish_clusters({KEY: weak}, 6, {}))


def test_the_boundary_is_inclusive():
    just_in = _cluster(1, 50.0, "in")
    just_out = _cluster(1, 49.9, "out")
    noise = [_cluster(3, 40.0, f"n{i}") for i in range(6)]
    selected = select_polish_clusters({KEY: noise + [just_in, just_out]}, 6, {}, strong_identity=50.0)
    assert (id(just_in), KEY) in selected
    assert (id(just_out), KEY) not in selected


def _cluster_of(key, n_genes, identity, contig):
    genes = ["sexP", "tptA", "rnhA", "glrA", "algA"][:n_genes]
    return types.SimpleNamespace(hits=[_hit(g, identity, contig, key=key) for g in genes])


def test_each_family_gets_its_own_tier_and_top_up():
    noise_a, true_a = _scenario()
    true_b = _cluster_of(OTHER, 1, 95.0, "b-true")
    noise_b = [_cluster_of(OTHER, 2, 38.0, f"b{i}") for i in range(7)]
    selected = select_polish_clusters({KEY: noise_a + [true_a], OTHER: noise_b + [true_b]}, 6, {},
                                      strong_identity=50.0)
    assert (id(true_a), KEY) in selected and (id(true_b), OTHER) in selected
    assert sum(1 for _, k in selected if k == KEY) == 6
    assert sum(1 for _, k in selected if k == OTHER) == 6


def test_strong_fragment_clusters_keep_their_priority_in_the_top_up():
    noise, true = _scenario()
    fragment = _cluster(1, 30.0, "frag")             # weak identity, but V3 strong
    strong = {(id(fragment), KEY): 200.0}
    selected = select_polish_clusters({KEY: noise + [true, fragment]}, 6, strong, strong_identity=50.0)
    assert (id(true), KEY) in selected               # identity tier
    assert (id(fragment), KEY) in selected           # first of the rest
    assert len(selected) == 6


# --- through the pipeline -------------------------------------------------

HITS = [
    ("mfa1", "c1", 100, 200, 40.0), ("pra1", "c1", 300, 400, 40.0), ("flk1", "c1", 500, 600, 40.0),
    ("mfa1", "c2", 100, 200, 90.0), ("pra1", "c2", 300, 400, 90.0),
]


def _run(tmp_path, **kw):
    _write_order(tmp_path, ORDER)
    _write_record(tmp_path)
    polished = []

    def localize(*a, **k):
        return [_tblastn(g, c, s, e, identity=i) for g, c, s, e, i in HITS]

    def polish(*, gene_name, window, **k):
        polished.append(window[0])
        start, end = next((s, e) for g, c, s, e, _ in HITS if g == gene_name and c == window[0])
        return _model(gene_name, window[0], start, end)

    run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=localize, polish_with_exonerate=polish, polish_with_miniprot=polish,
        evidence_diagnostics_path=tmp_path / "ed.jsonl", **kw,
    )
    return set(polished)


def test_the_pipeline_polishes_the_strong_cluster_by_default(tmp_path):
    # cap 1: c1 has more genes (3) at low identity, c2 fewer (2) at 90%.
    assert _run(tmp_path, max_polished_clusters_per_family=1) == {"c2"}


def test_turning_the_tier_off_restores_the_plain_cap(tmp_path):
    assert _run(tmp_path, max_polished_clusters_per_family=1, polish_strong_identity=None) == {"c1"}


def test_the_default_polishes_every_strong_cluster_even_past_the_cap(tmp_path):
    """The cap fixture of test_polish_cluster_cap: c1 3 genes at 80%, c2 2 genes at 70%, c3 2 genes at
    40%. With a cap of 1 the plain cap polishes c1 only; both c1 and c2 are at or above 50%, so the tier
    polishes both and leaves only the weak c3 capped."""
    from tests.detect.test_polish_cluster_cap import _run as run_cap_fixture
    (tmp_path / "plain").mkdir(); (tmp_path / "tiered").mkdir()
    _, plain = run_cap_fixture(tmp_path / "plain", max_polished_clusters_per_family=1, polish_strong_identity=None)
    _, tiered = run_cap_fixture(tmp_path / "tiered", max_polished_clusters_per_family=1)
    assert plain == {"c1"}
    assert tiered == {"c1", "c2"}


def test_with_the_cap_off_the_tier_changes_nothing(tmp_path):
    assert _run(tmp_path, max_polished_clusters_per_family=None) == {"c1", "c2"}


# --- command line -----------------------------------------------------------

def test_the_cli_default_and_zero():
    from MATPredict.__main__ import build_parser
    from MATPredict.detect.cli import _polish_strong_identity_from_args
    parser = build_parser()
    assert _polish_strong_identity_from_args(parser.parse_args(["detect", "--genome", "g.fa"])) == 50.0
    off = parser.parse_args(["detect", "--genome", "g.fa", "--polish-strong-identity", "0"])
    assert _polish_strong_identity_from_args(off) is None
    custom = parser.parse_args(["detect", "--genome", "g.fa", "--polish-strong-identity", "65"])
    assert _polish_strong_identity_from_args(custom) == 65.0


def test_an_args_object_without_the_flag_gets_the_default():
    from MATPredict.detect.cli import _polish_strong_identity_from_args
    assert _polish_strong_identity_from_args(types.SimpleNamespace()) == DEFAULT_POLISH_STRONG_IDENTITY
