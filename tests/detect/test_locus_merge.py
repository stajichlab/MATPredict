"""Calls from different families that describe the same physical locus are
reported once.

Curator's ruling 2026-09-27 (results/2026-09-27_caax_precursor/NOTE.md): the
Schizophyllum commune B locus was reported three times at one span -- once by
the generic Basidiomycota:PR family (receptor + CAAX precursor) and once each
by the curated Balpha and Bbeta families. One locus is reported that keeps
every family's evidence. Only families that share a roster `merge_group` may
merge, never HD with PR, and never calls on different contigs or with
conflicting idiomorph labels (e.g. the Cercospora kikuchii double call).
"""
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.locus_merge import merge_overlapping
from MATPredict.detect.pipeline import DetectionResult, GeneEvidence

PR = FamilyKey("Basidiomycota", "PR")
BA = FamilyKey("Basidiomycota", "Balpha")
BB = FamilyKey("Basidiomycota", "Bbeta")
HD = FamilyKey("Basidiomycota", "HD")
MAT = FamilyKey("Ascomycota", "MAT")
TUB = FamilyKey("Ascomycota", "MATtub")

GROUPS = {PR: ("B", True), BA: ("B", False), BB: ("B", False), HD: ("A", True)}


def _ev(gene, start, end, contig="c1"):
    return GeneEvidence(gene_name=gene, role="core_MAT", contig=contig, start=start,
                        end=end, strand="+", identity=90.0, coverage=None,
                        reference_record_id="r", method="exonerate_refine")


def _r(key, start=1000, end=9000, conf="high", idio="undetermined", contig="c1",
       genes=None):
    genes = genes or [("g", start, start + 100)]
    return DetectionResult(
        family_key=key, contig=contig, start=start, end=end, confidence=conf,
        idiomorph=idio, ambiguous_with=[], genes_found=[g for g, _, _ in genes],
        genes_missing=[], fragmented=False,
        gene_evidence=[_ev(g, s, e, contig) for g, s, e in genes],
    )


def test_the_schizophyllum_b_locus_is_reported_once():
    calls = [
        _r(PR, conf="medium", genes=[("caax_precursor", 2000, 2200),
                                     ("pheromone_receptor", 5000, 6000)]),
        _r(BA, genes=[("bar3", 5000, 6100), ("bap3-1", 6500, 6700)]),
        _r(BB, conf="medium", genes=[("bbr2", 3000, 4000), ("bbp2-1", 4100, 4300)]),
    ]
    out = merge_overlapping(calls, GROUPS)
    assert len(out) == 1
    m = out[0]
    assert m.family_key == BA  # specific family, higher confidence than Bbeta
    assert m.confidence == "high"
    fams = {f["family"] for f in m.merged_from}
    assert fams == {"Basidiomycota:PR", "Basidiomycota:Balpha", "Basidiomycota:Bbeta"}
    genes = {e.gene_name for e in m.gene_evidence}
    assert {"caax_precursor", "pheromone_receptor", "bar3", "bap3-1", "bbr2", "bbp2-1"} <= genes
    assert set(m.genes_found) >= {"bar3", "bbr2", "caax_precursor"}


def test_hd_never_merges_with_pr():
    out = merge_overlapping([_r(HD), _r(PR, conf="medium")], GROUPS)
    assert len(out) == 2
    assert all(r.merged_from == [] for r in out)


def test_different_contigs_never_merge():
    out = merge_overlapping([_r(PR), _r(BA, contig="c2")], GROUPS)
    assert len(out) == 2


def test_little_span_overlap_does_not_merge():
    # overlap 1,000 bp of a 10,000 bp shorter span = 10% < 50%
    out = merge_overlapping([_r(PR, 1000, 11000), _r(BA, 10000, 30000)], GROUPS)
    assert len(out) == 2


def test_conflicting_idiomorphs_do_not_merge():
    out = merge_overlapping([_r(BA, idio="3"), _r(BB, idio="2"),], {
        BA: ("B", False), BB: ("B", False)})
    assert {r.family_key for r in out} == {BA, BB}
    assert all(r.merged_from == [] for r in out)


def test_undetermined_is_compatible_and_the_determined_label_wins():
    out = merge_overlapping([_r(PR, conf="medium"), _r(BA, conf="medium", idio="3")], GROUPS)
    assert len(out) == 1
    assert out[0].idiomorph == "3"


def test_families_without_a_merge_group_are_untouched():
    # Ascomycota MAT / MATtub duplicates use different idiomorph vocabularies
    # and are not in a merge group (decision recorded in the roster note).
    a, b = _r(MAT, idio="MAT1-2"), _r(TUB, idio="MAT1-2")
    out = merge_overlapping([a, b], GROUPS)
    assert out == [a, b]


def test_the_merge_is_recorded_with_each_members_label_and_span():
    out = merge_overlapping([_r(PR, 1000, 9000, conf="medium"), _r(BA, 1500, 8000)], GROUPS)
    m = out[0]
    assert (m.start, m.end) == (1000, 9000)
    pr = next(f for f in m.merged_from if f["family"] == "Basidiomycota:PR")
    assert pr["confidence"] == "medium" and pr["start"] == 1000 and pr["end"] == 9000
    assert pr["idiomorph"] == "undetermined"


def test_order_of_unmerged_calls_is_kept():
    calls = [_r(HD), _r(PR, 50000, 60000, conf="medium"), _r(BA, 50500, 59000)]
    out = merge_overlapping(calls, GROUPS)
    assert [r.family_key for r in out] == [HD, BA]


def test_the_real_roster_groups_the_basidiomycota_b_and_a_loci_only():
    from pathlib import Path
    from MATPredict.detect.family_registry import load_all_families as load_families
    db = Path(__file__).resolve().parents[2] / "db"
    fams = {f.key: f for f in load_families(db)}
    groups = {k: (f.merge_group, f.merge_generic) for k, f in fams.items() if f.merge_group}
    assert groups[PR] == ("B", True)
    assert groups[BA] == ("B", False) and groups[BB] == ("B", False)
    assert groups[HD] == ("A", True)
    assert groups[FamilyKey("Basidiomycota", "Aalpha")] == ("A", False)
    assert not any(k.phylum == "Ascomycota" for k in groups)
