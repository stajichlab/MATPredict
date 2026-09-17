"""End-to-end integration check for the localize-then-polish flow (Task 11).

Exercises the real routing/family/order.yml data for the real, accepted
Mycosarcoma maydis a-locus record through the FULL genome-only Stage 0-3
flow: `search_localize` (Stage 1, tblastn) locates the family's core_MAT
genes, then `polish_with_exonerate`/`polish_with_miniprot` (Stage 2) refine
each one. Both search/polish stages are stubbed -- no live tblastn,
exonerate or miniprot binary is invoked -- but the routing, family loading,
record-family indexing, clustering, scoring, polish classification and
tiering all run for real against the real `db/`.

Companion to `test_integration_real_record.py`, which exercises the FAST
path (a supplied proteome) against the same real record. This file
exercises the GENOME-ONLY path (no proteome), which is the path that
actually calls `search_localize` and the polish wrappers.
"""
from __future__ import annotations
from pathlib import Path

import pytest

from MATPredict.detect.family_registry import load_all_families, route
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.polish import ExonSpan, PolishModel
from MATPredict.detect.reference_fasta import build_reference_fasta
from MATPredict.detect.search import SearchHit

DB_ROOT = Path("db")

# Confirmed real path (re-verified for this task): this record lives under
# "Ustilaginales" (an order), not "Ustilaginaceae" (a family) -- see
# `find db/Basidiomycota -iname "5270_*"` and test_integration_real_record.py's
# own note about this same path.
REAL_RECORD = DB_ROOT / "Basidiomycota" / "Ustilaginales" / "5270_521_aLocus_a1"

RECORD_ID = "5270_521_aLocus_a1"

# Real curated coordinates from REAL_RECORD/metadata.yaml (U37795.1, 1-based):
# mfa1 1319-1856 (-), pra1 3604-4914 (-). Used only as realistic tblastn/polish
# stand-in coordinates on a synthetic contig name -- no real genome FASTA is
# read by any stubbed call in this test.
GENE_COORDS = {"mfa1": (1319, 1856), "pra1": (3604, 4914)}


def _a_locus_family_key(families):
    return next(f.key for f in families if f.key.locus_name == "aLocus")


def _tblastn_hit(family_key, gene_name, start, end, identity=75.0):
    return SearchHit(
        family_key=family_key, gene_name=gene_name, role="core_MAT",
        contig="U37795.1", start=start, end=end, strand="-",
        identity=identity, reference_record_id=RECORD_ID, method="tblastn_genome",
    )


def _model(family_key, gene_name, start, end, *, method, identity=98.0):
    return PolishModel(
        gene_name=gene_name, family_key=family_key, role="core_MAT",
        contig="U37795.1", start=start, end=end, strand="-",
        exons=[ExonSpan(start, end)], identity=identity,
        reference_record_id=RECORD_ID, method=method,
    )


def _run(tmp_path, families, family_key, reference_fasta, *, miniprot_shift):
    """Genome-only run: `search_localize` finds both real genes, and both
    polish tools confirm them -- shifted by `miniprot_shift` bp to control
    whether the pair agrees (0) or disagrees (past the 10bp tolerance)."""

    def fake_localize(genome_fasta, localize_families, ref_fasta, record_families, runner=None):
        return [_tblastn_hit(family_key, name, *span) for name, span in GENE_COORDS.items()]

    def fake_exonerate(*, gene_name, **kwargs):
        start, end = GENE_COORDS[gene_name]
        return _model(family_key, gene_name, start, end, method="exonerate_refine")

    def fake_miniprot(*, gene_name, **kwargs):
        start, end = GENE_COORDS[gene_name]
        return _model(
            family_key, gene_name, start + miniprot_shift, end + miniprot_shift,
            method="miniprot_refine",
        )

    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa",
        proteome_fasta=None,
        taxid=5270,
        db_root=DB_ROOT,
        reference_fasta=reference_fasta,
        search_localize=fake_localize,
        polish_with_exonerate=fake_exonerate,
        polish_with_miniprot=fake_miniprot,
    )


@pytest.mark.skipif(not REAL_RECORD.exists(),
                     reason="requires the real curated record to be present")
def test_genome_only_pipeline_reaches_high_confidence_with_polished_agreement(tmp_path):
    families = route(5270, load_all_families(DB_ROOT))
    family_key = _a_locus_family_key(families)
    reference_fasta = build_reference_fasta(DB_ROOT, tmp_path / "reference.faa")

    outcome = _run(tmp_path, families, family_key, reference_fasta, miniprot_shift=0)
    a_locus_result = next(r for r in outcome.results if r.family_key.locus_name == "aLocus")

    assert a_locus_result.confidence == "high"
    assert "mfa1" in a_locus_result.genes_found
    assert "pra1" in a_locus_result.genes_found

    evidence = {e.gene_name: e for e in a_locus_result.gene_evidence}
    assert {"mfa1", "pra1"} <= set(evidence)
    # Both tools agree within tolerance -> status is the real computed
    # STATUS_AGREE value, never the GeneEvidence dataclass's placeholder default.
    assert evidence["mfa1"].status == "polished_agree"
    assert evidence["pra1"].status == "polished_agree"
    assert evidence["mfa1"].alternate_model is None
    assert a_locus_result.reference_records == [RECORD_ID]


@pytest.mark.skipif(not REAL_RECORD.exists(),
                     reason="requires the real curated record to be present")
def test_genome_only_pipeline_disagreeing_polish_models_still_reach_same_tier(tmp_path):
    """Integration-level version of test_pipeline.py's
    test_polished_agree_and_disagree_produce_identical_tier: run the real
    aLocus record's genes through the real routing/scoring/tiering stack
    twice, once with agreeing polish models and once with mfa1's miniprot
    model shifted 500bp away (past the 10bp tolerance -> polished_disagree),
    and confirm both reach the SAME confidence tier -- agreement is reported
    but never consulted by tiering."""
    families = route(5270, load_all_families(DB_ROOT))
    family_key = _a_locus_family_key(families)
    reference_fasta = build_reference_fasta(DB_ROOT, Path(str(tmp_path)) / "reference.faa")

    agree_outcome = _run(tmp_path / "agree", families, family_key, reference_fasta, miniprot_shift=0)
    disagree_outcome = _run(
        tmp_path / "disagree", families, family_key, reference_fasta, miniprot_shift=500
    )

    agree = next(r for r in agree_outcome.results if r.family_key.locus_name == "aLocus")
    disagree = next(r for r in disagree_outcome.results if r.family_key.locus_name == "aLocus")

    assert agree.confidence == disagree.confidence == "high"

    agree_evidence = {e.gene_name: e for e in agree.gene_evidence}
    disagree_evidence = {e.gene_name: e for e in disagree.gene_evidence}
    assert agree_evidence["mfa1"].status == "polished_agree"
    assert disagree_evidence["mfa1"].status == "polished_disagree"
    assert disagree_evidence["mfa1"].alternate_model is not None
    # canonical (exonerate) coordinates are identical either way -- exonerate
    # is preferred as canonical regardless of what miniprot reports.
    assert (agree_evidence["mfa1"].start, agree_evidence["mfa1"].end) == (
        disagree_evidence["mfa1"].start, disagree_evidence["mfa1"].end
    )
