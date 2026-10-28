"""scripts/idiomorph_labels.py and label_loci_table.py: the registry labels every Ascomycota and Mucoromycotina locus of the v0.7.0 tables; Basidiomycota stay pending; the grammar."""
import csv
import importlib.util
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]


def _load(name):
    sys.path.insert(0, str(_ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location(name, _ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


L = _load("idiomorph_labels")
T = _load("label_loci_table")


def _loci(campaign):
    with open(_ROOT / "results" / f"2026-10-26_v070_{campaign}" / "loci.tsv", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_grammar_examples():
    assert L.label_locus("Ascomycota:PM", "M")["label"] == "Taphr__SchPom_M"
    assert L.label_locus("Ascomycota:MTL", "A")["label"] == "Sacch__CanAlb_a"
    assert L.label_locus("Ascomycota:MATsc", "MATalpha")["label"] == "Sacch__SacCer_alpha"
    assert L.label_locus("Ascomycota:MATyl", "B")["label"] == "Sacch__YarLip_B"
    assert L.label_locus("Ascomycota:MATdc", "A")["label"] == "Sacch__GeoCan_A"
    assert L.label_locus("Mucoromycota:MAT", "Minus")["label"] == "Mucor__PhyBla_minus"
    assert L.label_locus("Ascomycota:MAT", "MAT1-2", "Eurotiomycetes")["label"] == "Peziz__AspFum_MAT1-2"
    assert L.label_locus("Ascomycota:MAT", "undetermined", "Arthoniomycetes")["label"] == "Peziz__NoRef_undet"


def test_gene_labels():
    assert L.gene_label("MTLalpha2", family="Ascomycota:MTL") == "CanAlb_alpha2"
    assert L.gene_label("MATa1", family="Ascomycota:MATog") == "OgaPol_a1"
    assert L.gene_label("MATa1") == "MATa1"            # KomPha and OgaPol: ambiguous without the family
    assert L.gene_label("matMc", "similarity_only") == "SchPom_Mc_like"
    assert L.gene_label("MAT1-1-1", family="Ascomycota:MAT", class_="Sordariomycetes") == "NeuCra_MAT1-1-1"
    assert L.gene_label("SLA2", family="Ascomycota:MAT") == "SLA2"      # flank gene


def test_every_ascomycota_and_mucoromycotina_locus_is_labelled():
    for campaign in ("ascomycota", "mucoromycotina"):
        rows = T.label_rows(_loci(campaign))
        assert {r["label_status"] for r in rows} == {"labelled"}, campaign


def test_basidiomycota_are_pending_and_keep_raw_names():
    rows = T.label_rows(_loci("basidiomycota"))
    assert {r["label_status"] for r in rows} == {"pending"}
    assert all(r["idiomorph_label"] == r["idiomorph"] for r in rows)


def test_registry_systems_are_declared_and_lineages_are_five_letters():
    reg = L.registry()
    assert all(len(k) == 5 for k in reg["lineages"])
    for fam, f in reg["families"].items():
        if f.get("system"): assert f["system"] in reg["systems"], fam
        assert f["lineage"] in reg["lineages"], fam
    assert set(reg["pezizomycotina_class_systems"].values()) <= set(reg["systems"])
