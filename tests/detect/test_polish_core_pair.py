"""Polish selection: a cluster with two core genes of one idiomorph is not capped out.

Seven Komagataella genomes have a real two-gene MTL cluster (MTLA1 + MTLA2, a MATa locus on a 2.1-2.3 kb
contig) that was never polished once nine families were searched: 18-24 MTL clusters were admitted
instead of 15-18, and the plain rank (distinct genes, then identity) put the small cluster outside the
cap of 6 (`polish_capped`). The tier makes the cluster's own property decide, not how many others there are.
"""
import json
import types

from MATPredict.detect.pipeline import (
    core_idiomorph_genes,
    has_core_pair,
    run_pipeline,
    select_polish_clusters,
)

from tests.detect.test_pipeline import _model, _tblastn, _write_order, _write_record

KEY = ("P", "MAT")
OTHER = ("P", "OTHER")
CORE = {KEY: {"mfa1": "a1", "pra1": "a1", "mfb1": "b1"}, OTHER: {"x1": "a1", "x2": "a1"}}


def _hit(gene, key=KEY, role="core_MAT", identity=40.0, superseded_by=None):
    return types.SimpleNamespace(gene_name=gene, identity=identity, family_key=key, role=role,
                                 superseded_by=superseded_by, method="tblastn_genome")


def _cluster(genes, key=KEY, identity=40.0):
    return types.SimpleNamespace(hits=[_hit(g, key, identity=identity) for g in genes])


def _noise(n):
    # 3 genes (more than the true cluster's 2), mixed idiomorphs, so no pair
    return [_cluster(["mfa1", "mfb1", "flk1"]) for _ in range(n)]


def test_two_core_genes_of_one_idiomorph_is_a_pair():
    assert has_core_pair(_cluster(["mfa1", "pra1"]), KEY, CORE)


def test_what_is_not_a_pair():
    assert not has_core_pair(_cluster(["mfa1"]), KEY, CORE)                       # one gene
    assert not has_core_pair(_cluster(["mfa1", "mfb1"]), KEY, CORE)               # two idiomorphs
    assert not has_core_pair(_cluster(["mfa1", "mfa1"]), KEY, CORE)               # one gene twice
    assert not has_core_pair(_cluster(["mfa1", "flk1"]), KEY, CORE)               # a flank is not a core gene
    assert not has_core_pair(_cluster(["x1", "x2"], key=OTHER), KEY, CORE)        # the other family's genes
    cluster = types.SimpleNamespace(hits=[_hit("mfa1"), _hit("pra1", superseded_by="mfb1")])
    assert not has_core_pair(cluster, KEY, CORE)                                  # a superseded hit is not a gene
    assert not has_core_pair(types.SimpleNamespace(hits=[_hit("mfa1"), _hit("pra1", role="flanking_conserved")]),
                             KEY, CORE)


def test_core_idiomorph_genes_reads_single_idiomorph_core_genes_only():
    family = types.SimpleNamespace(key=KEY, genes=[
        {"name": "A1", "role": "core_MAT", "present_in_idiomorphs": ["A"]},
        {"name": "both", "role": "core_MAT", "present_in_idiomorphs": ["A", "alpha"]},
        {"name": "none", "role": "core_MAT"},
        {"name": "flank", "role": "flanking_conserved", "present_in_idiomorphs": ["A"]},
    ])
    assert core_idiomorph_genes([family]) == {KEY: {"A1": "A"}}


def test_a_capped_out_pair_cluster_is_now_polished():
    true = _cluster(["mfa1", "pra1"])
    members = _noise(8) + [true]
    assert (id(true), KEY) not in select_polish_clusters({KEY: members}, 6, {}, None)   # the failure
    selected = select_polish_clusters({KEY: members}, 6, {}, None, CORE)
    assert (id(true), KEY) in selected
    assert len(selected) == 6                       # the pair takes the place of the weakest, same work


def test_unrelated_clusters_are_unchanged():
    members = _noise(9) + [_cluster(["mfa1"]), _cluster(["mfb1", "flk1"])]
    assert (select_polish_clusters({KEY: members}, 6, {}, 50.0, CORE)
            == select_polish_clusters({KEY: members}, 6, {}, 50.0))


def test_it_composes_with_the_identity_tier():
    strong = _cluster(["flk1"], identity=95.0)
    pair = _cluster(["mfa1", "pra1"])
    members = _noise(8) + [strong, pair]
    selected = select_polish_clusters({KEY: members}, 6, {}, 50.0, CORE)
    assert {(id(strong), KEY), (id(pair), KEY)} <= selected
    assert len(selected) == 6


def test_with_many_extra_families_the_selection_of_a_family_does_not_change():
    true = _cluster(["mfa1", "pra1"])
    members = _noise(8) + [true]
    alone = select_polish_clusters({KEY: members}, 6, {}, 50.0, CORE)
    crowded = {KEY: members, **{("P", f"F{i}"): _noise(20) for i in range(7)}, OTHER: _noise(12)}
    chosen = select_polish_clusters(crowded, 6, {}, 50.0, CORE)
    assert {s for s in chosen if s[1] == KEY} == alone
    assert (id(true), KEY) in chosen


def test_more_pair_clusters_than_the_cap_are_all_polished():
    pairs = [_cluster(["mfa1", "pra1"]) for _ in range(8)]
    selected = select_polish_clusters({KEY: pairs + _noise(4)}, 6, {}, None, CORE)
    assert selected == {(id(c), KEY) for c in pairs}


# --- through the pipeline ---------------------------------------------------

ORDER = (
    "phylum: P\nloci:\n"
    "  - locus_name: aLocus\n    vocabulary_type: pattern\n"
    "    idiomorph_pattern: \"^[ab][0-9]+$\"\n    taxonomic_scope: [1]\n"
    "    genes:\n"
    "      - {name: mfa1, role: core_MAT, present_in_idiomorphs: [a1]}\n"
    "      - {name: pra1, role: core_MAT, present_in_idiomorphs: [a1]}\n"
    "      - {name: mfb1, role: core_MAT, present_in_idiomorphs: [b1]}\n"
    "      - {name: flk1, role: flanking_conserved}\n"
)
HITS = ([(g, f"n{i}", 100 + 200 * j, 200 + 200 * j, 60.0)
         for i in range(3) for j, g in enumerate(["mfa1", "mfb1", "flk1"])]
        + [("mfa1", "true", 100, 200, 40.0), ("pra1", "true", 300, 400, 40.0)])


def _run(tmp_path, **kw):
    _write_order(tmp_path, ORDER)
    _write_record(tmp_path)
    polished = []

    def localize(*a, **k):
        return [_tblastn(g, c, s, e, identity=i) for g, c, s, e, i in HITS]

    def polish(*, gene_name, window, **k):
        polished.append(window[0])
        found = next(((s, e) for g, c, s, e, _ in HITS if g == gene_name and c == window[0]), None)
        if found is None:
            return None
        start, end = found
        return _model(gene_name, window[0], start, end)

    run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=localize, polish_with_exonerate=polish, polish_with_miniprot=polish,
        evidence_diagnostics_path=tmp_path / "ed.jsonl", **kw)
    rows = [json.loads(l) for l in (tmp_path / "ed.jsonl").read_text().splitlines()]
    return set(polished), {r["contig"]: r["polish_capped"] for r in rows if r["kind"] == "evidence"}


def test_the_pipeline_polishes_the_capped_out_pair_by_default(tmp_path):
    polished, capped = _run(tmp_path, max_polished_clusters_per_family=2, polish_strong_identity=None)
    assert "true" in polished and not capped["true"]
    assert len(polished) == 2                      # the pair takes a slot from the weakest ranked cluster


def test_the_flag_off_restores_the_earlier_cap(tmp_path):
    polished, capped = _run(tmp_path, max_polished_clusters_per_family=2, polish_strong_identity=None,
                            polish_core_pair=False)
    assert "true" not in polished and capped["true"]
    assert len(polished) == 2


def test_the_cli_flag():
    from MATPredict.__main__ import build_parser
    parser = build_parser()
    assert parser.parse_args(["detect", "--genome", "g.fa"]).polish_core_pair is True
    assert parser.parse_args(["detect", "--genome", "g.fa", "--no-polish-core-pair"]).polish_core_pair is False
