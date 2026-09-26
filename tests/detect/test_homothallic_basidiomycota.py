"""`homothallic_candidate` must not fire on partner genes of one locus.

Curator's ruling 2026-09-26: "the basidiomycetes agaricomycotina aren't really
homothallic in the same way so it doesn't really apply to try give that label".
The Basidiomycota anchor pilot (docs/notes/2026-09-26_basidiomycota-anchors.md,
branch basidio-anchors) measured the rule firing on 8 of 30 genomes: HD1+HD2
and bE+bW are partner subunits that sit together at one locus in EVERY strain,
and they carry different gene_classes, so the relaxed two-unrelated-genes path
took them for two idiomorphs.

Two changes, both data-driven:

1. A core gene with no `present_in_idiomorphs` is in every idiomorph, so it is
   evidence of no particular idiomorph and cannot make half of a homothallic
   pair. This is the general form of the HD1/HD2 case.
2. `homothallic_screen: false` in `order.yml` (phylum level, or per locus)
   turns the label off for a family. Basidiomycota sets it phylum-wide, which
   also covers Cryptococcus SXI1/SXI2, whose idiomorph restrictions are real.

The ascomycete behaviour is unchanged: D. hansenii and the CTG-clade a+alpha
loci, and Hydnotrya, still get the label.
"""
from pathlib import Path

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey, load_all_families
from MATPredict.detect.idiomorph import LOCUS_CLASS_HOMOTHALLIC, classify_locus
from MATPredict.detect.search import SearchHit

DB = Path(__file__).resolve().parents[2] / "db"


def _family(genes, screen=True, key=FamilyKey("Basidiomycota", "HD")):
    return Family(key, "pattern", None, "^A[0-9]+$", genes, [5338], homothallic_screen=screen)


HD = _family([
    {"name": "HD1", "role": "core_MAT", "gene_class": "HD1"},
    {"name": "HD2", "role": "core_MAT", "gene_class": "HD2"},
])
PR = _family([
    {"name": "STE3", "role": "core_MAT", "gene_class": "pheromone_receptor"},
    {"name": "phb", "role": "core_MAT", "gene_class": "pheromone_precursor"},
], key=FamilyKey("Basidiomycota", "PR"))
SXI = [
    {"name": "SXI1", "role": "core_MAT", "gene_class": "HD1", "present_in_idiomorphs": ["alpha"]},
    {"name": "SXI2", "role": "core_MAT", "gene_class": "HD2", "present_in_idiomorphs": ["a"]},
]


def _hit(family, gene, start, end, method="tblastn_genome"):
    return SearchHit(family.key, gene, "core_MAT", "c1", start, end, "+", 60.0, "rec1", method)


def _cluster(*hits):
    return GeneCluster("c1", min(h.start for h in hits), max(h.end for h in hits), list(hits))


def test_hd1_and_hd2_at_one_locus_are_not_homothallic():
    """Both full-length, different classes, 300 bp apart: the pilot's case."""
    cluster = _cluster(_hit(HD, "HD1", 1_000, 2_800), _hit(HD, "HD2", 3_100, 4_900))
    got = classify_locus(cluster, HD, full_length_models=frozenset({"HD1", "HD2"}))
    assert got != LOCUS_CLASS_HOMOTHALLIC


def test_partner_genes_found_in_the_proteome_are_not_homothallic_either():
    cluster = _cluster(_hit(HD, "HD1", 1_000, 2_800, "diamond_proteome"),
                       _hit(HD, "HD2", 3_100, 4_900, "diamond_proteome"))
    assert classify_locus(cluster, HD) != LOCUS_CLASS_HOMOTHALLIC


def test_a_receptor_and_a_pheromone_are_not_homothallic():
    cluster = _cluster(_hit(PR, "STE3", 1_000, 2_200), _hit(PR, "phb", 3_000, 3_200))
    got = classify_locus(cluster, PR, full_length_models=frozenset({"STE3", "phb"}))
    assert got != LOCUS_CLASS_HOMOTHALLIC


def test_the_screen_switch_turns_the_label_off():
    """SXI1/SXI2 are idiomorph-restricted: only the switch stops them."""
    on, off = _family(SXI, screen=True), _family(SXI, screen=False)
    cluster = _cluster(_hit(on, "SXI1", 1_000, 2_800), _hit(on, "SXI2", 3_100, 4_900))
    models = frozenset({"SXI1", "SXI2"})
    assert classify_locus(cluster, on, full_length_models=models) == LOCUS_CLASS_HOMOTHALLIC
    assert classify_locus(cluster, off, full_length_models=models) != LOCUS_CLASS_HOMOTHALLIC


def test_the_curated_database_turns_it_off_for_basidiomycota_only():
    families = load_all_families(DB)
    basidio = [f for f in families if f.key.phylum == "Basidiomycota"]
    assert basidio and not any(f.homothallic_screen for f in basidio)
    assert all(f.homothallic_screen for f in families if f.key.phylum != "Basidiomycota")


def test_a_locus_setting_overrides_the_phylum_setting(tmp_path):
    (tmp_path / "P").mkdir()
    (tmp_path / "P" / "order.yml").write_text(
        "phylum: P\nhomothallic_screen: false\nloci:\n"
        "  - {locus_name: A, vocabulary_type: enum, idiomorph_values: [x, y],"
        " taxonomic_scope: [1], genes: [], homothallic_screen: true}\n"
        "  - {locus_name: B, vocabulary_type: enum, idiomorph_values: [x, y],"
        " taxonomic_scope: [1], genes: []}\n"
    )
    screens = {f.key.locus_name: f.homothallic_screen for f in load_all_families(tmp_path)}
    assert screens == {"A": True, "B": False}


def test_debaryomyces_hansenii_is_still_a_homothallic_candidate():
    """Regression, CBS767: MTLA2 (HMG box, A) and MTLalpha1 (alpha box, alpha),
    both full-length models ~1 kb apart, on the curated Ascomycota:MTL family."""
    [mtl] = [f for f in load_all_families(DB) if f.key == FamilyKey("Ascomycota", "MTL")]
    cluster = _cluster(_hit(mtl, "MTLA2", 10_000, 10_900), _hit(mtl, "MTLalpha1", 11_900, 12_600))
    got = classify_locus(cluster, mtl, full_length_models=frozenset({"MTLA2", "MTLalpha1"}))
    assert got == LOCUS_CLASS_HOMOTHALLIC


def test_hydnotrya_is_still_a_homothallic_candidate():
    [mat] = [f for f in load_all_families(DB) if f.key == FamilyKey("Ascomycota", "MATtub")]
    cluster = _cluster(_hit(mat, "MAT1-1-1", 8_931, 9_867), _hit(mat, "MAT1-2-1", 11_840, 12_418))
    got = classify_locus(cluster, mat, full_length_models=frozenset({"MAT1-1-1", "MAT1-2-1"}))
    assert got == LOCUS_CLASS_HOMOTHALLIC
