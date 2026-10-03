"""V3 polish-cap rank: strong-fragment clusters first, inside the cap.

Curator's ruling 2026-09-30 (results/2026-09-30_cap_protection/NOTE.md): within
the per-family polish cap, a cluster whose unmodelled tblastn HSP fragment
already scores at or above the family classifier's gate threshold -- and is not
claimed by a paralog class -- is ranked FIRST. The cap size never grows.
Families without a classifier keep the plain rank. Measured by replay: V3
recovered 19 of 20 displaced calls in a stress test (caps 3 and 2 against a
cap-off reference) with no wrong cluster promoted, and recovered M. pusillus
NRRL A-13674 at cap 6.
"""
import types

from MATPredict.detect.classifier import ClassifierVerdict
from MATPredict.detect.pipeline import (
    _polish_rank,
    rank_for_polish,
    select_polish_clusters,
    strong_fragment_scores,
)

KEY = ("P", "MAT")


def _hit(gene, identity, contig="c1", key=KEY, method="tblastn_genome"):
    return types.SimpleNamespace(gene_name=gene, identity=identity, family_key=key,
                                 method=method, contig=contig, start=1, end=90, strand="+")


def _cluster(n_genes, identity, contig):
    genes = ["sexP", "tptA", "rnhA", "glrA", "algA"][:n_genes]
    return types.SimpleNamespace(hits=[_hit(g, identity, contig) for g in genes])


def _members():
    # Seven clusters: six 3-gene noise clusters, and one 1-gene true locus that
    # the plain rank puts last.
    noise = [_cluster(3, 35.0, f"n{i}") for i in range(6)]
    true = _cluster(1, 30.0, "true")
    return noise, true


def test_a_strong_fragment_cluster_is_ranked_first_within_the_cap():
    noise, true = _members()
    members = noise + [true]
    strong = {(id(true), KEY): 224.8}
    plain = sorted(members, key=lambda c: _polish_rank(c, KEY))
    assert plain[-1] is true  # the plain rank leaves it outside a cap of 6
    ranked = rank_for_polish(members, KEY, strong)
    assert ranked[0] is true
    selected = select_polish_clusters({KEY: members}, 6, strong)
    assert (id(true), KEY) in selected
    assert len(selected) == 6


def test_the_cap_size_is_never_exceeded_even_when_many_clusters_are_strong():
    clusters = [_cluster(2, 40.0, f"s{i}") for i in range(9)]
    strong = {(id(c), KEY): 150.0 for c in clusters}
    selected = select_polish_clusters({KEY: clusters}, 6, strong)
    assert len(selected) == 6


def test_without_strong_clusters_the_selection_is_the_plain_rank():
    noise, true = _members()
    members = noise + [true]
    selected = select_polish_clusters({KEY: members}, 6, {})
    plain = {(id(c), KEY) for c in sorted(members, key=lambda c: _polish_rank(c, KEY))[:6]}
    assert selected == plain


def _family(key=KEY, classifier=True):
    genes = [{"name": "sexP", "role": "core_MAT", "present_in_idiomorphs": ["Plus"]},
             {"name": "sexM", "role": "core_MAT", "present_in_idiomorphs": ["Minus"]},
             {"name": "tptA", "role": "flanking_conserved"}]
    spec = {"type": "hmm", "dir": "/nonexistent", "min_margin": 25} if classifier else None
    return types.SimpleNamespace(key=key, genes=genes, idiomorph_classifier=spec)


def _verdict(scores, paralog_class=None):
    return ClassifierVerdict(scores=scores, margin=abs(scores["Plus"] - scores["Minus"]),
                             idiomorph="Plus", min_margin=25.0, proteins_scored=1,
                             classifier_input="hsp_fragment", paralog_class=paralog_class)


def _strong(members, family, verdict):
    return strong_fragment_scores(
        {family.key: members}, [family], genome_fasta=None, genetic_code=1,
        load_fn=lambda spec, fam: object(), gate_fn=lambda spec: (100.0, "manifest"),
        translate_fn=lambda *a, **k: "MKTAYIAKQR", classify_fn=lambda *a, **k: verdict,
    )


def test_a_fragment_at_or_above_the_gate_is_strong():
    _, true = _members()
    out = _strong([true], _family(), _verdict({"Plus": 224.8, "Minus": 40.0}))
    assert out == {(id(true), KEY): 224.8}


def test_a_fragment_below_the_gate_is_not_strong():
    _, true = _members()
    out = _strong([true], _family(), _verdict({"Plus": 86.4, "Minus": 30.0}))
    assert out == {}


def test_a_paralog_claimed_fragment_is_not_promoted():
    _, true = _members()
    out = _strong([true], _family(), _verdict({"Plus": 30.0, "Minus": 120.0}, paralog_class="P1"))
    assert out == {}


def test_a_family_without_a_classifier_gets_no_strong_clusters():
    _, true = _members()
    out = _strong([true], _family(classifier=False), _verdict({"Plus": 224.8, "Minus": 40.0}))
    assert out == {}


def test_only_the_familys_own_core_gene_hsps_are_scored():
    seen = []
    fam = _family()
    c = types.SimpleNamespace(hits=[
        _hit("tptA", 90.0),                      # flank: not a classifier gene
        _hit("sexP", 90.0, key=("Q", "MAT")),    # another family's hit
        _hit("sexP", 90.0, method="miniprot"),   # not a tblastn HSP
    ])
    strong_fragment_scores(
        {fam.key: [c]}, [fam], genome_fasta=None, genetic_code=1,
        load_fn=lambda spec, f: object(), gate_fn=lambda spec: (100.0, "manifest"),
        translate_fn=lambda *a, **k: seen.append(1) or "MK",
        classify_fn=lambda *a, **k: _verdict({"Plus": 500.0, "Minus": 0.0}),
    )
    assert seen == []
