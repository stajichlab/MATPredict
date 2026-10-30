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
        assert {r["label_mapping"] for r in rows} == {"labelled"}, campaign


def test_basidiomycota_labels():
    assert L.label_locus("Basidiomycota:HD", "undetermined")["label"] == "Basid__CopCin_HD"
    assert L.label_locus("Basidiomycota:HD", "A43")["label"] == "Basid__CopCin_HD-A43"
    assert L.label_locus("Basidiomycota:redPR", "A1")["label"] == "Basid__RhoTor_PR-A1"
    assert L.label_locus("Basidiomycota:wallMAT", "v2")["label"] == "Basid__WalMel_MAT-v2"
    assert L.label_locus("Basidiomycota:bLocus", "undetermined")["label"] == "Basid__UstMay_b"
    assert L.label_locus("Basidiomycota:MAT", "alpha")["label"] == "Basid__CryDen_alpha"
    assert L.gene_label("bE", family="Basidiomycota:rustHD") == "PucGra_bE"
    assert L.gene_label("MIP1", family="Basidiomycota:HD") == "MIP1"


def test_every_basidiomycota_locus_is_labelled():
    rows = T.label_rows(_loci("basidiomycota"))
    assert {r["label_mapping"] for r in rows} == {"labelled"}


def test_registry_systems_are_declared_and_lineages_are_five_letters():
    reg = L.registry()
    assert all(len(k) == 5 for k in reg["lineages"])
    for fam, f in reg["families"].items():
        if f.get("system"): assert f["system"] in reg["systems"], fam
        assert f["lineage"] in reg["lineages"], fam
    assert set(reg["pezizomycotina_class_systems"].values()) <= set(reg["systems"])


# ---- registry version and frozen codes (#210, curator 2026-10-29)
import copy
import csv as _csv
import json

import yaml


def test_registry_version_is_date_plus_hash_and_the_stamp_is_current():
    v = L.registry_version()
    date, h = v.split("+")
    assert len(date) == 10 and len(h) == 8 and h == L.registry_hash()
    assert L.registry_problems() == []            # a rule edited without `--stamp` fails here


def test_registry_hash_changes_with_a_rule_and_not_with_the_comment():
    reg = copy.deepcopy(L.registry())
    base = L.registry_hash(reg)
    assert L.registry_hash(copy.deepcopy(reg)) == base                      # same registry, same version
    reg["_comment"] += " edited"
    assert L.registry_hash(reg) == base
    reg["families"]["Ascomycota:MATsc"]["cells"]["MATa"] = "A"            # a rule change
    assert L.registry_hash(reg) != base
    assert any("registry_hash" in p for p in L.registry_problems(reg))


def test_frozen_codes_cannot_be_removed_or_renamed_without_an_alias():
    reg = copy.deepcopy(L.registry())
    assert {"Peziz", "Sacch", "Taphr", "Mucor", "Basid"} <= set(reg["frozen_codes"]["lineages"])
    assert {"SacCer", "CanAlb", "KomPha", "OgaPol", "PicKud", "YarLip", "GeoCan", "SchPom", "PhyBla", "CopCin", "SchCom", "CryDen", "UstMay", "RhoTor", "PucGra", "WalMel",
            "NeuCra", "AspFum", "BipMay", "BotCin", "XanPar", "TubMel", "NoRef"} <= set(reg["frozen_codes"]["systems"])
    # removing a frozen code fails
    gone = copy.deepcopy(reg); del gone["systems"]["OgaPol"]
    assert any("OgaPol" in p and "alias" in p for p in L.registry_problems(gone))
    # a rename needs an alias, and the alias target has to exist
    renamed = copy.deepcopy(reg)
    renamed["systems"]["OgaNew"] = renamed["systems"].pop("OgaPol")
    assert any("OgaPol" in p for p in L.registry_problems(renamed))
    renamed["aliases"]["systems"]["OgaPol"] = "OgaNew"
    assert not [p for p in L.registry_problems(renamed) if "OgaPol" in p]
    renamed["aliases"]["systems"]["OgaPol"] = "Missing"
    assert any("Missing" in p for p in L.registry_problems(renamed))
    # an alias resolves to the new code
    L.registry()["aliases"]["systems"]["__t"] = "SacCer"
    try:
        assert L.canonical_code("systems", "__t") == "SacCer" and L.canonical_code("systems", "SacCer") == "SacCer"
    finally:
        del L.registry()["aliases"]["systems"]["__t"]


def test_every_curated_family_has_a_registry_entry():
    reg = L.registry()["families"]
    for f in sorted((_ROOT / "db").glob("*/order.yml")):
        for locus in yaml.safe_load(f.read_text())["loci"]:
            assert f"{f.parent.name}:{locus['locus_name']}" in reg, (f.parent.name, locus["locus_name"])
    assert L.label_locus("Ascomycota:MATpk", "MATalpha")["label"] == "Sacch__PicKud_alpha"
    assert L.label_locus("Ascomycota:mat3", "M")["label"] == "Taphr__SchPom_M"
    assert L.label_locus("Basidiomycota:Abeta", "undetermined")["label"] == "Basid__SchCom_Abeta"


# ---- label_status (#211)
def test_status_precedence_and_reasons():
    st = L.locus_status("Ascomycota:MATsc", "undetermined", cross_call=True)
    assert st == {"label_status": "cross_call", "label_reasons": "cross_call;undetermined"}
    st = L.locus_status("Ascomycota:MTL", "A", family_unconfirmed=True, cross_call=True)
    assert st == {"label_status": "family_unconfirmed", "label_reasons": "family_unconfirmed;cross_call"}
    assert L.locus_status("Ascomycota:MAT", "MAT1-1", "Arthoniomycetes")["label_status"] == "no_reference"
    assert L.locus_status("Ascomycota:MATsc", "MATa") == {"label_status": "ok", "label_reasons": ""}
    assert L.locus_status("Ascomycota:MATsc", "undetermined")["label_status"] == "undetermined"
    assert L.locus_status("Basidiomycota:HD", "undetermined")["label_status"] == "ok"      # a locus label without an allele is complete
    assert L.locus_status("Ascomycota:MATsc", "MATa", similarity_only=True)["label_status"] == "similarity_only"
    assert L.STATUS_ORDER == ("family_unconfirmed", "cross_call", "no_reference", "similarity_only", "undetermined", "ok")


def _labelled(campaign):
    rows = _loci(campaign)
    raw = [(r["idiomorph"], r["family_called"]) for r in rows]
    T.label_rows(rows)
    assert [(r["idiomorph"], r["family_called"]) for r in rows] == raw      # no raw name or call changes
    return rows


def test_status_counts_on_the_v070_loci():
    asco, basid, muc = _labelled("ascomycota"), _labelled("basidiomycota"), _labelled("mucoromycotina")
    for rows in (asco, basid, muc):
        assert {r["label_status"] for r in rows} <= set(L.STATUS_ORDER)         # exactly one status per locus
        assert all(r["registry_version"] == L.registry_version() for r in rows)
    reasons = lambda rows, k: sum(k in r["label_reasons"].split(";") for r in rows)
    assert reasons(asco, "no_reference") == 63                                    # the NoRef loci of v0.7.0
    assert reasons(asco, "undetermined") + reasons(basid, "undetermined") + reasons(muc, "undetermined") == 932      # the `_undet` labels (#211: about 930)
    yeast = {r["genome"] for r in asco if {"Ascomycota:MATsc", "Ascomycota:MTL"} <= {x["family_called"] for x in asco if x["genome"] == r["genome"]}}
    assert len({r["genome"] for r in asco if "cross_call" in r["label_reasons"] and r["family_called"] in T.YEAST_SYSTEMS}) == 181
    assert yeast                                                                   # MATsc + MTL genomes are among them


def test_every_duplicate_family_call_of_the_default_dedupe_is_marked_cross_call():
    dd = _ROOT / "results" / "2026-10-28_idiomorph_crosscall_flags" / "loci_v070_ascomycota_deduped.tsv"
    dup = {(r["genome"], r["contig"], r["start"], r["end"], r["family_called"]) for r in _csv.DictReader(open(dd, newline=""), delimiter="\t") if r.get("duplicate_family_call") == "true"}
    assert len(dup) == 185
    rows = _labelled("ascomycota")
    got = {(r["genome"], r["contig"], r["start"], r["end"], r["family_called"]) for r in rows if "cross_call" in r["label_reasons"]}
    assert dup <= got
