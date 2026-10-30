"""The MATpk (Pichia kudriavzevii) family and its two CBS 573 records, curator decision 2026-10-29 (#189)."""
from pathlib import Path

import yaml

from MATPredict.detect.family_registry import FamilyKey, load_all_families, route

DB = Path(__file__).resolve().parents[2] / "db"
PICHIA_GENUS, PICHIA_KUDRIAVZEVII = 4919, 4909
BASE = [131567, 2759, 33154, 4751, 451864, 4890, 716545, 147537, 3239874, 3243775, 1156497]


def _family():
    return next(f for f in load_all_families(DB) if f.key == FamilyKey("Ascomycota", "MATpk"))


def _record(rid):
    return yaml.safe_load((DB / "Ascomycota" / "Pichiales" / rid / "metadata.yaml").read_text())


def test_matpk_declaration():
    fam = _family()
    assert list(fam.taxonomic_scope) == [PICHIA_KUDRIAVZEVII]   # the species, not the genus (replay, README)
    assert fam.fallback_searchable is False and fam.fallback_subphylum == 147537
    # heterothallic: the second allele is not a cassette, and no homothallism claim is made
    assert fam.cassette_report is False and fam.two_idiomorphs_report is False
    assert set(fam.idiomorph_classes.values()) == {"unassigned"}
    genes = {g["name"]: g for g in fam.genes}
    assert set(genes) == {"MATa1", "MATa2", "MATalpha1", "MATalpha2", "sla2", "tgl1"}
    assert genes["MATalpha1"]["gene_class"] == "alpha_box" and genes["MATa2"].get("gene_class") is None
    assert genes["sla2"].get("optional") and genes["tgl1"].get("optional")


def test_matpk_records_are_one_strain_two_alleles():
    a, m = _record("4909_cbs-573_MATpk_MATa"), _record("4909_cbs-573_MATpk_MATalpha")
    for r in (a, m):
        assert r["mating_type"]["locus_name"] == "MATpk" and r["mating_type"]["system"] == "heterothallic"
        assert r["taxonomy"]["taxid"] == 4909 and r["organism"]["strain"]["name"] == "CBS 573"
        assert r["validation"]["status"] == "accepted" and r["validation"]["taxonomy_current"] is True
        assert {g["status"] for g in r["validation"]["sequence_match"]["per_gene"]} == {"pass"}
    assert a["mating_type"]["idiomorphs"] == ["MATa"] and m["mating_type"]["idiomorphs"] == ["MATalpha"]
    assert a["locus"]["assembly_accession"] == "GCA_003054445.1" and m["locus"]["assembly_accession"] is None


def test_matpk_mata2_is_the_two_exon_model_with_the_62bp_intron():
    genes = {g["name"]: g for g in _record("4909_cbs-573_MATpk_MATa")["genes"]}
    mata2 = genes["MATa2"]
    assert mata2["strand"] == "+" and mata2["exons"] == [{"start": 311575, "end": 311841}, {"start": 311904, "end": 312245}]
    assert 311904 - 311841 - 1 == 62          # the intron the RNA-seq splice reads support
    proteins = (DB / "Ascomycota" / "Pichiales" / "4909_cbs-573_MATpk_MATa" / "proteins.faa").read_text().split(">")[1:]
    mata2_protein = next(p for p in proteins if "name=MATa2" in p.splitlines()[0])
    assert len("".join(mata2_protein.splitlines()[1:])) == 202
    assert genes["MATa1"]["strand"] == "-" and (genes["MATa1"]["start"], genes["MATa1"]["end"]) == (310491, 311300)


def test_pichia_kudriavzevii_routes_to_matpk_and_other_pichia_species_do_not():
    fams = load_all_families(DB)
    d = route(4909, fams, lineage_taxids_resolver=lambda t: BASE + [PICHIA_GENUS, 4909], phylum_name_resolver=lambda t: "Ascomycota")
    assert d.routing_mode == "direct" and {f.key.locus_name for f in d.families} == {"MATpk"}
    # a strain under the species routes by lineage too
    d = route(1234567, fams, lineage_taxids_resolver=lambda t: BASE + [PICHIA_GENUS, 4909, 1234567], phylum_name_resolver=lambda t: "Ascomycota")
    assert d.routing_mode == "lineage" and {f.key.locus_name for f in d.families} == {"MATpk"}
    # another Pichia species (P. kluyveri 36015) falls back to the phylum and does not search MATpk (scope-only, never in fallback)
    d = route(36015, fams, lineage_taxids_resolver=lambda t: BASE + [PICHIA_GENUS, 36015], phylum_name_resolver=lambda t: "Ascomycota")
    assert d.routing_mode == "phylum_fallback" and "MATpk" not in {f.key.locus_name for f in d.families}
    # another Pichiaceae genus is not routed here and does not search it in the phylum fallback
    d = route(36023, fams, lineage_taxids_resolver=lambda t: BASE + [461281, 36023], phylum_name_resolver=lambda t: "Ascomycota")
    assert "MATpk" not in {f.key.locus_name for f in d.families}
