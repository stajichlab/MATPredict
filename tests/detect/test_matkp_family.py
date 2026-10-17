"""The MATkp (Komagataella) family declaration, curator decisions of 2026-10-17/18.

gene_class is declared on MATalpha1 (alpha_box) and MATa2 (HMG_box) only; MATa1 and MATalpha2 are homeodomain
proteins and the schema has no homeodomain class. idiomorph_class stays unassigned and two_idiomorphs_report
stays off (open curator questions).
"""
from pathlib import Path

from MATPredict.detect.family_registry import FamilyKey, load_all_families

DB = Path(__file__).resolve().parents[2] / "db"


def _family():
    return next(f for f in load_all_families(DB) if f.key == FamilyKey("Ascomycota", "MATkp"))


def test_matkp_gene_classes():
    classes = {g["name"]: g.get("gene_class") for g in _family().genes}
    assert classes["MATa2"] == "HMG_box"
    assert classes["MATalpha1"] == "alpha_box"
    assert classes["MATa1"] is None
    assert classes["MATalpha2"] is None


def test_matkp_open_questions_stay_open():
    fam = _family()
    assert set(fam.idiomorph_classes.values()) == {"unassigned"}
    assert fam.two_idiomorphs_report is False
    assert fam.fallback_searchable is False


# ---- the Saccharomycotina MATa2 profile is supporting evidence only (curator, 2026-10-18)
import shutil
import subprocess
import sys

import pytest
import yaml

PROFILE = DB / "Ascomycota" / "classifiers" / "MATa2_saccharomycotina"
REPO = Path(__file__).resolve().parents[2]


def test_mata2_profile_is_not_referenced_anywhere():
    # order.yml COMMENTS may mention the profile as evidence; no declared value may point at it.
    for phylum_order in DB.glob("*/order.yml"):
        doc = yaml.safe_load(phylum_order.read_text())
        assert "MATa2_saccharomycotina" not in yaml.safe_dump(doc), phylum_order
        for locus in doc["loci"]:
            assert "idiomorph_classifier" not in locus or "MATa2" not in str(locus["idiomorph_classifier"])
    for py in (REPO / "src").rglob("*.py"):
        assert "MATa2_saccharomycotina" not in py.read_text(), py


def test_mata2_manifest_states_its_limits():
    m = yaml.safe_load((PROFILE / "manifest.yaml").read_text())
    text = m["status_statement"]
    assert "NOT referenced" in text and "incomplete" in text and "thin" in text
    assert len(m["training"]) == 7
    assert all(t["source_accession"] and t["species"] for t in m["training"])


@pytest.mark.skipif(shutil.which("mafft") is None, reason="mafft not on PATH")
def test_mata2_profile_builds_reproducibly(tmp_path):
    outs = []
    for k in range(2):
        d = tmp_path / f"b{k}"
        d.mkdir()
        shutil.copy(PROFILE / "training.faa", d / "training.faa")
        subprocess.run([sys.executable, str(REPO / "scripts" / "build_mata2_saccharomycotina_hmm.py"),
                        "--out", str(d)], check=True, env={**__import__("os").environ, "PYTHONPATH": str(REPO / "src")})
        outs.append(yaml.safe_load((d / "manifest.yaml").read_text()))
    assert outs[0]["hmm_sha256"] == outs[1]["hmm_sha256"]
    committed = yaml.safe_load((PROFILE / "manifest.yaml").read_text())
    assert outs[0]["training_sha256"] == committed["training_sha256"]
    if outs[0]["tool_versions"] == committed["tool_versions"]:
        assert outs[0]["hmm_sha256"] == committed["hmm_sha256"]


# ---- Ogataea lives in its own scope-only family MATog (curator, 2026-10-18)
def _og():
    return next(f for f in load_all_families(DB) if f.key == FamilyKey("Ascomycota", "MATog"))


def test_matkp_scope_stays_komagataella_and_holds_no_ogataea_record():
    assert list(_family().taxonomic_scope) == [460517]  # NOT widened to Ogataea 461281 or Pichiaceae 1156497
    ids = {p.name for p in (DB / "Ascomycota" / "Pichiales").iterdir() if p.is_dir()}
    assert not any("MATkp" in i and i.startswith(("460523_", "1005962_")) for i in ids)


def test_matog_family_declaration():
    fam = _og()
    assert list(fam.taxonomic_scope) == [461281]  # genus Ogataea (taxonkit-verified)
    assert fam.fallback_searchable is False
    assert fam.two_idiomorphs_report is False
    assert set(fam.idiomorph_classes.values()) == {"unassigned"}
    genes = {g["name"]: g for g in fam.genes}
    assert genes["MATa2"].get("optional") is True and genes["MATa2"].get("gene_class") is None
    assert genes["MATalpha1"]["gene_class"] == "alpha_box"
    assert genes["MATa1"].get("gene_class") is None and genes["MATalpha2"].get("gene_class") is None
    assert set(genes) == {"MATa1", "MATa2", "MATalpha1", "MATalpha2", "sla2", "dic1"}


def test_matog_records():
    ids = {"460523_cbs-4732_MATog_MATa-MATalpha", "1005962_dl-1_MATog_MATa-MATalpha"}
    for rid in ids:
        rec = yaml.safe_load((DB / "Ascomycota" / "Pichiales" / rid / "metadata.yaml").read_text())
        assert rec["mating_type"]["locus_name"] == "MATog"
        assert rec["mating_type"]["system"] == "pseudohomothallic"
        assert rec["validation"]["status"] == "accepted"
        assert "MATa2" not in {g["name"] for g in rec["genes"]}  # not recorded as present or absent (unresolved)
