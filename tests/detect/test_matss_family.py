"""The MATss (Saccharomycopsis schoenii) family and its two S. schoenii candidate records (curator, 2026-10-30)."""
from pathlib import Path

import yaml

from MATPredict.detect.family_registry import FamilyKey, load_all_families, route

DB = Path(__file__).resolve().parents[2] / "db"


def _family():
    return next(f for f in load_all_families(DB) if f.key == FamilyKey("Ascomycota", "MATss"))


def test_matss_declaration():
    fam = _family()
    assert list(fam.taxonomic_scope) == [44080]
    assert fam.fallback_searchable is False and fam.cassette_report is False and fam.two_idiomorphs_report is False
    genes = {g["name"]: g for g in fam.genes}
    assert genes["MATalpha1"]["gene_class"] == "alpha_box"
    assert genes["sla2"]["gene_class"] == "sla2_homolog"
    assert set(genes) == {"MATa1", "MATa2", "MATalpha1", "MATalpha2", "sla2", "dic1"}


def test_matss_records_are_accepted_tier_3_models_with_table_12():
    for rid in ("44080_gdd12_MATss_MATa-MATalpha", "44080_gdd12_MATss_MATa"):
        rec = yaml.safe_load((DB / "Ascomycota" / "Ascoideales" / rid / "metadata.yaml").read_text())
        assert rec["mating_type"]["locus_name"] == "MATss" and rec["mating_type"]["system"] == "homothallic"
        assert rec["evidence"]["boundaries"]["tier"] == 3
        assert all(g["transl_table"] == 12 for g in rec["genes"])
        assert rec["validation"]["status"] == "accepted"
        assert (DB / "Ascomycota" / "Ascoideales" / rid / "proteins.faa").stat().st_size > 0


def test_matss_routes_the_species_only():
    fams = load_all_families(DB)
    base = [131567, 2759, 33154, 4751, 451864, 4890, 716545, 147537, 4891, 2926619, 34366, 4943]
    d = route(44080, fams, lineage_taxids_resolver=lambda t: base + [44080], phylum_name_resolver=lambda t: "Ascomycota")
    assert d.routing_mode in ("direct", "lineage") and "MATss" in {f.key.locus_name for f in d.families}
    # another Saccharomycopsis species (S. fibuligera 4944, not in scope) is not routed to MATss
    d = route(4944, fams, lineage_taxids_resolver=lambda t: base + [4944], phylum_name_resolver=lambda t: "Ascomycota")
    assert "MATss" not in {f.key.locus_name for f in d.families}
    d = route(4909, fams, lineage_taxids_resolver=lambda t: [131567, 2759, 33154, 4751, 451864, 4890, 716545, 147537, 4909],
              phylum_name_resolver=lambda t: "Ascomycota")
    assert "MATss" not in {f.key.locus_name for f in d.families}


def test_matbz_and_matpf_records_are_accepted_tier_3_models():
    ids = {"Phaffomycetales": ["36038_pb4207_MATbz_MATa-MATalpha", "36022_joy008_MATpf_MATa", "36022_tbrc-4498_MATpf_MATalpha"]}
    for order, rids in ids.items():
        for rid in rids:
            rec = yaml.safe_load((DB / "Ascomycota" / order / rid / "metadata.yaml").read_text())
            assert rec["mating_type"]["locus_name"] in ("MATbz", "MATpf")
            assert rec["evidence"]["boundaries"]["tier"] == 3
            assert rec["validation"]["status"] == "accepted"
            assert (DB / "Ascomycota" / order / rid / "proteins.faa").stat().st_size > 0
