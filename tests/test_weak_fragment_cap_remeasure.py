"""Fixture test of scripts/weak_fragment_cap_remeasure.py (#341)."""
from __future__ import annotations

import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("wfc_remeasure", Path(__file__).resolve().parents[1] / "scripts" / "weak_fragment_cap_remeasure.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def ev(gene, status, bits):
    return {"gene": gene, "status": status, "bitscore": bits}


def locus(idiom, tier, genes):
    return {"family": "Ascomycota:MAT", "contig": "c", "start": 1, "end": 9, "idiomorph": idiom, "confidence": tier, "gene_evidence": genes}


def test_old_caps_a_call_with_a_thin_other_idiomorph_fragment_and_new_does_not(tmp_path):
    docs = [
        ("g1", {"detected": [locus("MAT1-1", "high", [ev("MAT1-1-1", "polished_agree", 300.0), ev("MAT1-2-1", "unpolished", 40.0)])]}),
        ("g2", {"detected": [locus("MAT1-1", "medium", [ev("MAT1-1-1", "unpolished", 41.0)])]}),
        ("g3", {"detected": [locus("MAT1-2", "low", [ev("MAT1-2-1", "unpolished", 41.0)])]}),
        ("g4", {"detected": [{"family": "Ascomycota:MATsc", "idiomorph": "MAT1-1", "confidence": "high", "gene_evidence": []}]}),
    ]
    # g1: MAT1-1-1 agrees and is the strongest gene under both rules, so neither caps; g2: both cap; g3 is already low.
    c = mod.run(docs, tmp_path)
    rows = (tmp_path / "changes.tsv").read_text().splitlines()
    assert [r.split("\t")[0] for r in rows[1:]] == ["g2"]
    assert c["old"]["matched"] == 1 and c["new"]["matched"] == 1
    assert c["old"]["loci"] == 3 and c["new"]["loci"] == 3          # MATsc is out of scope


def test_the_rules_differ_when_the_other_idiomorph_gene_is_polished():
    # MAT1-2-1 polished_single at 40 bits outranks an unpolished MAT1-1-1 at 300 bits under the old rule only.
    old, new = mod.evaluate_both(locus("MAT1-1", "high", [ev("MAT1-1-1", "unpolished", 300.0), ev("MAT1-2-1", "polished_single", 40.0)]))
    assert old is not None and old["strongest_core_gene"] == "MAT1-2-1"
    assert new is None


def test_a_directory_of_reports_is_read(tmp_path):
    import yaml
    d = tmp_path / "runs" / "g2"
    d.mkdir(parents=True)
    (d / "detection_report.yaml").write_text(yaml.safe_dump({"detected": [locus("MAT1-1", "medium", [ev("MAT1-1-1", "unpolished", 41.0)])]}))
    assert mod.main([str(tmp_path / "runs"), str(tmp_path / "out")]) == 0
    assert (tmp_path / "out" / "changes.tsv").read_text().splitlines()[1].startswith("g2\t")
