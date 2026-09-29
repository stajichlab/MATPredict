"""A noise-level hit that could not be modelled must not cap a call at medium.

Review finding F6 (results/2026-09-28_fable_review/README.md), curator's
ruling 2026-09-28. gzUmbRama1: an rnhA HSP at 38.3% / 34.7 bits, 28.4 kb from
the locus, fell inside the 50 kb cluster gap; neither tool modelled it, and
`any_gene_unpolished` capped a 100%-identity, 4-modelled-gene locus at medium
(the record itself says rnhA is not at the Umbelopsis locus). Same in Umbra1.

Rule: a gene whose every hit in the call scores below the family's noise floor
(`flank_carried_min_bitscore`, 39 bits -- the measured known-noise maximum was
30.4 bits, results/2026-09-27_flank_bitscore_floor/) is not a gene "that
failed to model"; it is ignored by the unpolished cap, whatever its role. A
gene with any hit at or above the floor, or with no bitscore recorded (e.g. a
proteome fast-path hit), still counts.
"""
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import sub_floor_genes
from MATPredict.detect.search import SearchHit

KEY = FamilyKey("Mucoromycota", "MAT")
OTHER = FamilyKey("Ascomycota", "MAT")


def _hit(gene, bits, key=KEY, role="flanking_conserved"):
    return SearchHit(key, gene, role, "c1", 1, 100, "+", 38.0, "rec1",
                     "tblastn_genome", bitscore=bits)


def test_a_gene_whose_every_hit_is_below_the_floor_is_sub_floor():
    assert sub_floor_genes([_hit("rnhA", 34.7)], KEY, 39.0) == {"rnhA"}


def test_one_hit_at_or_above_the_floor_keeps_the_gene():
    assert sub_floor_genes([_hit("rnhA", 34.7), _hit("rnhA", 39.0)], KEY, 39.0) == set()


def test_a_hit_without_a_bitscore_keeps_the_gene():
    assert sub_floor_genes([_hit("rnhA", None)], KEY, 39.0) == set()


def test_core_genes_are_covered_too():
    assert sub_floor_genes([_hit("sexM", 25.0, role="core_MAT")], KEY, 39.0) == {"sexM"}


def test_other_families_hits_are_ignored():
    assert sub_floor_genes([_hit("rnhA", 20.0, key=OTHER)], KEY, 39.0) == set()
