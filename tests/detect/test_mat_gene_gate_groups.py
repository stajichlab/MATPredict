"""Taxon-scoped flank support in the MAT-gene gate (curator's ruling 2026-10-01).

In Circinella, Thamnostylum, Fennellomyces and Zychaea the sex gene sits 2-9 kb
from rnhA and tptA/algA/glrA sit elsewhere, so a real locus can never show the
default two roster flanks. The roster names a `flank_support_groups` entry: for
a genome whose taxid lineage meets the group's taxids, rnhA alone is flank
support. Every other genome -- other Mucorales, Lichtheimia, a run without a
taxid -- keeps the default rule, unchanged.
"""
from pathlib import Path

import yaml

from MATPredict.detect.family_registry import FamilyKey, load_all_families
from MATPredict.detect.mat_gene_gate import (
    WITHHELD_MAT_GENE_GATE, apply_mat_gene_gate, flank_support_group,
)
from MATPredict.detect.pipeline import DetectionResult, GeneEvidence
from MATPredict.detect.report import _result_doc

KEY = FamilyKey("Mucoromycota", "MAT")
GROUP = {"name": "circinella_group", "taxids": [101102, 64650, 101135, 64651], "flanks": ["rnhA"]}
SPEC = {"type": "hmm", "dir": "/x", "min_margin": 25, "flank_support_groups": [GROUP]}
SPEC_NO_GROUP = {"type": "hmm", "dir": "/x", "min_margin": 25}

#: Lineages (taxid + ancestors). 4827 = Mucorales, 499202 = Lichtheimiaceae.
FENNELLOMYCES = frozenset([1329386, 101135, 499202, 4827])
LICHTHEIMIA = frozenset([688394, 688353, 499202, 4827])
RHIZOPUS = frozenset([64495, 4842, 4827])

DB = Path(__file__).resolve().parents[2] / "db"


def _ev(gene, role, identity, status="polished_single"):
    return GeneEvidence(gene, role, "c1", 1, 100, "+", identity, None, "rec1",
                        "exonerate_refine", status=status)


def _result(evidence, scores, idiomorph="Plus"):
    clf = {"method": "hmm", "idiomorph": idiomorph, "scores": scores,
           "margin": abs(scores["Plus"] - scores["Minus"]),
           "classifier_input": "model", "min_margin": 25}
    return DetectionResult(
        family_key=KEY, contig="c1", start=1, end=1000, confidence="medium",
        idiomorph=idiomorph, ambiguous_with=[],
        genes_found=sorted({e.gene_name for e in evidence}), genes_missing=[],
        fragmented=False, gene_evidence=evidence, polished_genes=2,
        idiomorph_classifier=clf,
    )


# Fennellomyces RSA 1415 shape: sexP 98.9 bits (below the gate), rnhA alone.
SEXP = _ev("sexP", "core_MAT", 45.0)
RNHA = _ev("rnhA", "flanking_conserved", 78.0)
FENN_CALL = _result([SEXP, RNHA], {"Plus": 98.9, "Minus": 30.3})


def test_rnhA_alone_keeps_a_call_in_the_group():
    kept, withheld = apply_mat_gene_gate([FENN_CALL], {KEY: SPEC}, FENNELLOMYCES)
    assert withheld == [] and len(kept) == 1
    assert kept[0].mat_gene_gate_group == {
        "group": "circinella_group", "flanks": ["rnhA"], "best_score": 98.9,
        "mat_gene_min_score": 100.0,
    }
    assert _result_doc(kept[0])["mat_gene_gate_group"]["group"] == "circinella_group"


def test_the_same_call_outside_the_group_is_withheld():
    for lineage in (RHIZOPUS, LICHTHEIMIA):
        kept, withheld = apply_mat_gene_gate([FENN_CALL], {KEY: SPEC}, lineage)
        assert kept == []
        assert [w.withheld_reason for w in withheld] == [WITHHELD_MAT_GENE_GATE]
        # No group applied, so the withheld detail is exactly the pre-group one.
        assert "flank_support_group" not in withheld[0].withheld_detail


def test_a_run_without_a_taxid_keeps_the_default_rule():
    kept, withheld = apply_mat_gene_gate([FENN_CALL], {KEY: SPEC}, None)
    assert kept == [] and len(withheld) == 1


def test_the_group_needs_its_own_flank_not_any_flank():
    r = _result([SEXP, _ev("glrA", "flanking_variable", 77.0)], {"Plus": 98.9, "Minus": 30.3})
    kept, withheld = apply_mat_gene_gate([r], {KEY: SPEC}, FENNELLOMYCES)
    assert kept == []
    assert withheld[0].withheld_detail["flank_support_group"] == "circinella_group"
    assert withheld[0].withheld_detail["flank_support_group_flanks"] == ["rnhA"]


def test_a_weak_or_unmodelled_rnhA_does_not_count():
    for ev in (_ev("rnhA", "flanking_conserved", 35.0),
               _ev("rnhA", "flanking_conserved", 80.0, status="unpolished")):
        r = _result([SEXP, ev], {"Plus": 98.9, "Minus": 30.3})
        kept, _ = apply_mat_gene_gate([r], {KEY: SPEC}, FENNELLOMYCES)
        assert kept == []


def test_calls_the_default_rules_keep_are_not_marked():
    r = _result([SEXP, RNHA], {"Plus": 164.3, "Minus": 30.7})
    kept, _ = apply_mat_gene_gate([r], {KEY: SPEC}, FENNELLOMYCES)
    assert kept == [r] and kept[0].mat_gene_gate_group is None
    assert "mat_gene_gate_group" not in _result_doc(kept[0])


def test_a_family_without_groups_ignores_the_lineage():
    kept, withheld = apply_mat_gene_gate([FENN_CALL], {KEY: SPEC_NO_GROUP}, FENNELLOMYCES)
    assert kept == [] and len(withheld) == 1


def test_a_paralog_class_is_still_withheld_first():
    r = _result([SEXP, RNHA], {"Plus": 98.9, "Minus": 30.3})
    r.idiomorph_classifier["paralog_class"] = "P1"
    r.idiomorph_classifier["paralog_scores"] = {"P1": 150.0}
    kept, withheld = apply_mat_gene_gate([r], {KEY: SPEC}, FENNELLOMYCES)
    assert kept == [] and withheld[0].withheld_reason == "paralog_class"


def test_the_shipped_roster_scopes_the_group_to_the_four_genera_only():
    family = next(f for f in load_all_families(DB) if f.key == KEY)
    groups = family.idiomorph_classifier["flank_support_groups"]
    assert [g["name"] for g in groups] == ["circinella_group"]
    assert sorted(groups[0]["taxids"]) == sorted([101102, 64650, 101135, 64651])
    assert groups[0]["flanks"] == ["rnhA"]
    spec = family.idiomorph_classifier
    # Lichtheimia (688353), the family (499202), Rhizomucor (4838) and the
    # order (4827) are not in the group.
    for taxid in (688353, 499202, 4838, 4827):
        assert flank_support_group(spec, frozenset([taxid])) is None
    for taxid in (101102, 64650, 101135, 64651):
        assert flank_support_group(spec, frozenset([taxid]))["name"] == "circinella_group"


def test_the_group_neighbours_are_not_roster_genes():
    """VMA1, pntA and gpmI are recorded in the two records' extended_flank
    only. A roster gene, even unsearched, enters other Mucorales' scoring
    denominators where `searchable_genes` is not applied."""
    doc = yaml.safe_load((DB / "Mucoromycota" / "order.yml").read_text())
    names = {g["name"] for g in doc["loci"][0]["genes"]}
    assert not names & {"VMA1", "pntA", "gpmI"}
