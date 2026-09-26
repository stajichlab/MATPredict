"""Two unrelated core genes, both modelled, at one locus are a homothallic candidate.

Curator's ruling, 2026-09-25 ("yes relax the homothallic rule").

The proteome-only guard exists because sexM and sexP are BOTH HMG-box genes:
on a genome-only run, two partial tblastn alignments of one HMG region looked
like two genes (54 of 75 false candidates in the 44-genus sweep). It stays for
such same-domain pairs. But it also blocked every genome-only run from ever
emitting the class, including pairs that cannot be one region seen twice:
MAT1-1-1 (alpha box) with MAT1-2-1 (HMG box), or MTLalpha1 (alpha box) with
MTLa2 (HMG box). Measured 2026-09-25: two Hydnotrya genomes with both
Pezizomycotina genes modelled 2-3 kb apart, and six CTG-clade genomes with
MTLa and MTLalpha1 genes at one locus, all labelled a single idiomorph.

The relaxed path needs all of: both genes FULL-LENGTH MODELS (passed in as
`full_length_models`; see `pipeline.full_length_models`), both carrying a
`gene_class`, the two classes DIFFERENT, and the two hits not overlapping.

"Full-length" and "not in a cross-match" are what make it safe. Simulated
first on the existing reports: a first version without them would have
labelled 113 loci homothallic in heterothallic-rich panels (Ophiostomatales
42/130, Diaporthales 25/216, Lecanorales 11, Dothideomycetes 11, ...). Two
real-biology causes: MAT1-1-3 (HMG box, in the MAT1-1 idiomorph) cross-matching
MAT1-2-1, and truncated MAT1-1-1 remnants inside MAT1-2 idiomorphs
(Ophiostomatales, e.g. a 266 bp fragment in Leptographium). Requiring each
model to cover >= 50% of its own reference protein and to have taken no part
in an idiomorph resolution removed all 113, and kept both Hydnotrya loci and
the CTG-clade a+alpha loci.
"""
from MATPredict.detect.idiomorph import IdiomorphResolution
from MATPredict.detect.pipeline import full_length_models
from MATPredict.detect.polish import (
    STATUS_AGREE, STATUS_UNPOLISHED, ExonSpan, PolishModel, PolishOutcome,
)
from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.idiomorph import LOCUS_CLASS_HOMOTHALLIC, classify_locus
from MATPredict.detect.search import SearchHit

PEZ = Family(
    FamilyKey("Ascomycota", "MATtub"), "enum", ["MAT1-1", "MAT1-2"], None,
    [
        {"name": "MAT1-1-1", "role": "core_MAT", "gene_class": "alpha_box",
         "present_in_idiomorphs": ["MAT1-1"]},
        {"name": "MAT1-2-1", "role": "core_MAT", "gene_class": "HMG_box",
         "present_in_idiomorphs": ["MAT1-2"]},
        {"name": "SLA2", "role": "flanking_conserved", "gene_class": "sla2_homolog"},
    ],
    [147549],
)
MUC = Family(
    FamilyKey("Mucoromycota", "MAT"), "enum", ["Plus", "Minus"], None,
    [
        {"name": "sexP", "role": "core_MAT", "present_in_idiomorphs": ["Plus"]},
        {"name": "sexM", "role": "core_MAT", "present_in_idiomorphs": ["Minus"]},
    ],
    [4827],
)


def _tbn(family, gene, start, end):
    return SearchHit(family.key, gene, "core_MAT", "c1", start, end, "+", 49.0,
                     "rec1", "tblastn_genome")


def _cluster(*hits):
    return GeneCluster("c1", min(h.start for h in hits), max(h.end for h in hits), list(hits))


def test_hydnotrya_both_genes_polished_is_homothallic():
    """GCA_040803455.1: MAT1-1-1 8,931-9,867 and MAT1-2-1 11,840-12,418."""
    cluster = _cluster(_tbn(PEZ, "MAT1-1-1", 8931, 9867), _tbn(PEZ, "MAT1-2-1", 11840, 12418))
    got = classify_locus(cluster, PEZ, full_length_models=frozenset({"MAT1-1-1", "MAT1-2-1"}))
    assert got == LOCUS_CLASS_HOMOTHALLIC


def test_without_the_models_it_stays_what_it_was():
    """The guard still holds when nothing was modelled."""
    cluster = _cluster(_tbn(PEZ, "MAT1-1-1", 8931, 9867), _tbn(PEZ, "MAT1-2-1", 11840, 12418))
    assert classify_locus(cluster, PEZ) != LOCUS_CLASS_HOMOTHALLIC


def test_only_one_gene_modelled_is_not_enough():
    cluster = _cluster(_tbn(PEZ, "MAT1-1-1", 8931, 9867), _tbn(PEZ, "MAT1-2-1", 11840, 12418))
    got = classify_locus(cluster, PEZ, full_length_models=frozenset({"MAT1-1-1"}))
    assert got != LOCUS_CLASS_HOMOTHALLIC


def test_overlapping_hits_are_not_two_genes():
    cluster = _cluster(_tbn(PEZ, "MAT1-1-1", 1000, 1600), _tbn(PEZ, "MAT1-2-1", 1500, 2100))
    got = classify_locus(cluster, PEZ, full_length_models=frozenset({"MAT1-1-1", "MAT1-2-1"}))
    assert got != LOCUS_CLASS_HOMOTHALLIC


def test_same_domain_pairs_keep_the_proteome_guard():
    """sexM/sexP carry no differing gene_class: polished models do not open
    the relaxed path, exactly the case the guard was built for."""
    cluster = _cluster(_tbn(MUC, "sexP", 1000, 1600), _tbn(MUC, "sexM", 4000, 4600))
    got = classify_locus(cluster, MUC, full_length_models=frozenset({"sexP", "sexM"}))
    assert got != LOCUS_CLASS_HOMOTHALLIC


def test_too_far_apart_is_still_not_one_locus():
    cluster = _cluster(_tbn(PEZ, "MAT1-1-1", 1000, 1600), _tbn(PEZ, "MAT1-2-1", 900000, 900600))
    got = classify_locus(cluster, PEZ, full_length_models=frozenset({"MAT1-1-1", "MAT1-2-1"}))
    assert got != LOCUS_CLASS_HOMOTHALLIC


def _outcome(gene, aa, status=STATUS_AGREE, record="rec1"):
    model = PolishModel(gene, PEZ.key, "core_MAT", "c1", 1, aa * 3, "+",
                        [ExonSpan(1, aa * 3)], 50.0, record, "exonerate_refine")
    return PolishOutcome(status=status, canonical=model, exonerate_model=model, miniprot_model=None)


def test_a_model_covering_half_its_reference_is_full_length():
    polish_by = {(7, PEZ.key, "MAT1-1-1"): _outcome("MAT1-1-1", 200)}
    got = full_length_models(polish_by, {7}, PEZ.key, {("rec1", "MAT1-1-1"): 400}, set())
    assert got == frozenset({"MAT1-1-1"})


def test_a_truncated_remnant_is_not_full_length():
    """Leptographium's MAT1-1-1 remnant inside a MAT1-2 idiomorph: 266 bp."""
    polish_by = {(7, PEZ.key, "MAT1-1-1"): _outcome("MAT1-1-1", 88)}
    got = full_length_models(polish_by, {7}, PEZ.key, {("rec1", "MAT1-1-1"): 400}, set())
    assert got == frozenset()


def test_a_gene_in_a_cross_match_resolution_is_excluded():
    """MAT1-2-1 that beat MAT1-1-3 at one position may BE MAT1-1-3."""
    polish_by = {(7, PEZ.key, "MAT1-2-1"): _outcome("MAT1-2-1", 300)}
    got = full_length_models(polish_by, {7}, PEZ.key, {("rec1", "MAT1-2-1"): 394},
                             {"MAT1-2-1", "MAT1-1-3"})
    assert got == frozenset()


def test_an_unpolished_gene_is_not_a_model():
    polish_by = {(7, PEZ.key, "MAT1-1-1"): _outcome("MAT1-1-1", 400, status=STATUS_UNPOLISHED)}
    got = full_length_models(polish_by, {7}, PEZ.key, {("rec1", "MAT1-1-1"): 400}, set())
    assert got == frozenset()


def test_another_clusters_model_does_not_count():
    polish_by = {(8, PEZ.key, "MAT1-1-1"): _outcome("MAT1-1-1", 400)}
    got = full_length_models(polish_by, {7}, PEZ.key, {("rec1", "MAT1-1-1"): 400}, set())
    assert got == frozenset()
