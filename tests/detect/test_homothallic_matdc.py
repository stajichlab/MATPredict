"""MATdc: MATA (HMG-like) and MATB (alpha-box-like) together at one locus are a
homothallic candidate (curator request 2026-10-17; the pattern of Geotrichum
ghanense GCA_022023575.1). Uses the real order.yml family, so the gene classes
declared there are what is tested."""
from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import FamilyKey, load_all_families
from MATPredict.detect.idiomorph import LOCUS_CLASS_HOMOTHALLIC, LOCUS_CLASS_MAT, classify_locus
from MATPredict.detect.search import SearchHit


def _family(real_db_root):
    return next(f for f in load_all_families(real_db_root) if f.key == FamilyKey("Ascomycota", "MATdc"))


def _hit(family, gene, role, start, end, strand="+"):
    return SearchHit(family.key, gene, role, "JAKKCX010000006.1", start, end, strand, 49.0,
                     "rec", "tblastn_genome")


def _cluster(*hits):
    return GeneCluster("c1", min(h.start for h in hits), max(h.end for h in hits), list(hits))


def test_ghanense_pattern_both_models_is_homothallic_candidate(real_db_root):
    fam = _family(real_db_root)
    # G. ghanense: MATB 1,653,133-1,653,648; MATA 1,667,948-1,669,624; sla2 1,672,203-1,675,406
    cluster = _cluster(
        _hit(fam, "MATB", "core_MAT", 1653133, 1653648),
        _hit(fam, "MATA", "core_MAT", 1667948, 1669624, "-"),
        _hit(fam, "sla2", "flanking_conserved", 1672203, 1675406, "-"),
    )
    got = classify_locus(cluster, fam, full_length_models=frozenset({"MATA", "MATB"}))
    assert got == LOCUS_CLASS_HOMOTHALLIC


def test_unmodelled_or_one_sided_calls_are_not_homothallic(real_db_root):
    fam = _family(real_db_root)
    both = _cluster(
        _hit(fam, "MATB", "core_MAT", 1653133, 1653648),
        _hit(fam, "MATA", "core_MAT", 1667948, 1669624, "-"),
        _hit(fam, "sla2", "flanking_conserved", 1672203, 1675406, "-"),
    )
    # tblastn-only hits with no full-length models: the guard holds
    assert classify_locus(both, fam) == LOCUS_CLASS_MAT
    # one model is full length, the other is not
    assert classify_locus(both, fam, full_length_models=frozenset({"MATB"})) == LOCUS_CLASS_MAT
    # a single idiomorph (the ordinary Geotrichum case) stays a mat_locus
    one = _cluster(_hit(fam, "MATA", "core_MAT", 584221, 585066),
                   _hit(fam, "sla2", "flanking_conserved", 579055, 582249))
    assert classify_locus(one, fam, full_length_models=frozenset({"MATA"})) == LOCUS_CLASS_MAT


def test_too_far_apart_is_not_one_locus(real_db_root):
    fam = _family(real_db_root)
    far = _cluster(_hit(fam, "MATB", "core_MAT", 1000, 1500),
                   _hit(fam, "MATA", "core_MAT", 1000 + 60000, 61900, "-"))
    assert classify_locus(far, fam, full_length_models=frozenset({"MATA", "MATB"})) != LOCUS_CLASS_HOMOTHALLIC
