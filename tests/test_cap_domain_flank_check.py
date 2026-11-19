"""Fixture test of scripts/cap_domain_flank_check.py (#341): flank distances, domain expectation and the verdict table."""
from __future__ import annotations

import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("cdf", Path(__file__).resolve().parents[1] / "scripts" / "cap_domain_flank_check.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def ev(gene, start, end, contig="c1", **kw):
    return {"gene": gene, "start": start, "end": end, "contig": contig, **kw}


def test_expected_domains():
    assert mod.expected_domain("MAT1-1-1") == "alpha-box"
    assert mod.expected_domain("MAT1-2-1") == "HMG-box" and mod.expected_domain("MAT1-1-3") == "HMG-box"
    assert mod.expected_domain("MAT1-1-5") is None and mod.expected_domain("MAT1-2-7") is None


def test_flank_offsets_are_signed_and_sorted_and_ignore_other_contigs_and_lost_genes():
    gene = ev("MAT1-2-1", 10000, 11000)
    locus = {"gene_evidence": [gene, ev("APN2", 2000, 3000), ev("SLA2", 15000, 16000), ev("COX13", 1, 2, contig="c2"), ev("SLA2", 50, 60, lost_to="x")]}
    assert mod.flank_offsets(locus, gene) == [("SLA2", 4000), ("APN2", -7000)]


def test_verdicts_without_genome_sequence_and_for_genes_with_no_expected_domain():
    gene = ev("MAT1-1-5", 10000, 11000)
    locus = {"gene_evidence": [gene, ev("SLA2", 15000, 16000)]}
    r = mod.check_locus(locus, {"strongest_own_gene": "MAT1-1-5"}, lambda c: None, [])
    assert r["domain"] == "n/a" and r["flank_ok"] and r["verdict"] == "flank_only"
    far = {"gene_evidence": [gene, ev("SLA2", 90000, 91000)]}
    assert mod.check_locus(far, {"strongest_own_gene": "MAT1-1-5"}, lambda c: None, [])["verdict"] == "neither"
    hmg = {"gene_evidence": [ev("MAT1-2-1", 10000, 11000), ev("SLA2", 15000, 16000)]}
    r = mod.check_locus(hmg, {"strongest_own_gene": "MAT1-2-1"}, lambda c: None, [])
    assert r["domain"] == "not_checked" and r["verdict"] == "not_checked"


def test_domain_found_in_a_real_orf(monkeypatch):
    # a stub of the domain search: one HMG-box ORF overlapping the model
    monkeypatch.setattr(mod, "domain_hits", lambda seq, off, hmms, min_aa: [{"acc": "PF00505", "domain": "HMG-box", "score": 40.0, "evalue": 1e-9, "orf_start": 10100, "orf_end": 10400, "aa": 100}])
    gene = ev("MAT1-2-1", 10000, 11000)
    locus = {"gene_evidence": [gene, ev("APN2", 2000, 3000)]}
    r = mod.check_locus(locus, {"strongest_own_gene": "MAT1-2-1"}, lambda c: "A" * 20000, [])
    assert r["domain"] == "found" and r["flank_ok"] and r["verdict"] == "domain+flank"
    far = {"gene_evidence": [gene]}
    assert mod.check_locus(far, {"strongest_own_gene": "MAT1-2-1"}, lambda c: "A" * 20000, [])["verdict"] == "domain_only"
    monkeypatch.setattr(mod, "domain_hits", lambda *a: [])
    assert mod.check_locus(locus, {"strongest_own_gene": "MAT1-2-1"}, lambda c: "A" * 20000, [])["verdict"] == "flank_only"
