"""A roster gene with its own `taxonomic_scope` (`--scoped-roster-policy`).

Prototype 2026-10-22: GAP1 and BGL1 sit 2.6-3.6 kb from the MAT gene in
Pleosporaceae (Bipolaris, Curvularia, Exserohilum), where SLA2/APN2 are
180-290 kb away, so the MAT cluster held one gene and the evidence floor
dropped it (analysis/2026-10-18_cochliobolus-mat-recall.md). Added to the
whole Ascomycota MAT roster they lowered 24 of 94 Dothideomycetes controls
from high to medium. Scoped, they reach only Pleosporaceae genomes, and only
with the policy on; everywhere else the routed roster and the query set are
the ones a run had before the genes were added.
"""
from pathlib import Path

import pytest

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import (
    Family, FamilyKey, load_all_families, load_record_families, route, scope_roster_genes,
)
from MATPredict.detect.reference_fasta import build_reference_fasta, searchable_genes_by_family
from MATPredict.detect.scoring import score_cluster
from MATPredict.detect.search import SearchHit

PLEOSPORACEAE, BIPOLARIS, BMAYDIS = 28556, 33194, 5016
DOTHIDEALES, AUREOBASIDIUM = 5195, 5579
LEOTIOMYCETA, DOTHIDEOMYCETES = 716546, 147541
SCOPED = {"GAP1", "BGL1"}


def _lineage(taxid):
    """Ancestors as NCBI gives them (the genome's own taxid is not included)."""
    if taxid == BMAYDIS:
        return [4890, 147538, LEOTIOMYCETA, DOTHIDEOMYCETES, 92860, PLEOSPORACEAE, BIPOLARIS]
    if taxid == AUREOBASIDIUM:
        return [4890, 147538, LEOTIOMYCETA, DOTHIDEOMYCETES, DOTHIDEALES]
    raise LookupError(taxid)


def _db_mat():
    return next(f for f in load_all_families(Path("db"))
                if f.key == FamilyKey("Ascomycota", "MAT"))


def _names(family):
    return [g["name"] for g in family.genes]


def _routed_mat(taxid, policy, **kw):
    d = route(taxid, load_all_families(Path("db")), lineage_taxids_resolver=_lineage,
              scoped_roster_policy=policy, **kw)
    return d, next(f for f in d.families if f.key == FamilyKey("Ascomycota", "MAT"))


def test_the_db_roster_declares_gap1_bgl1_pleosporaceae_only_as_optional_slots():
    genes = {g["name"]: g for g in _db_mat().genes}
    for name in SCOPED:
        g = genes[name]
        assert g["taxonomic_scope"] == [PLEOSPORACEAE]
        assert g["role"] == "flanking_conserved" and g["optional"] is True and g["slot"] == name


def test_policy_off_removes_the_scoped_genes_even_inside_the_scope():
    _, mat = _routed_mat(BMAYDIS, "off")
    assert not SCOPED & set(_names(mat))
    assert _names(mat) == [n for n in _names(_db_mat()) if n not in SCOPED]


def test_policy_on_keeps_them_inside_the_scope_only():
    _, inside = _routed_mat(BMAYDIS, "on")
    assert SCOPED <= set(_names(inside))
    _, outside = _routed_mat(AUREOBASIDIUM, "on")
    assert _names(outside) == [n for n in _names(_db_mat()) if n not in SCOPED]
    assert "GAP1" not in outside.gene_aliases and "GAP1" in inside.gene_aliases


def test_operator_routes_and_failed_lookups_remove_them():
    _, mat = _routed_mat(None, "on", phylum="Ascomycota")
    assert not SCOPED & set(_names(mat))
    d = route(4242, load_all_families(Path("db")), lineage_taxids_resolver=_lineage,
              phylum_name_resolver=lambda t: "Ascomycota", scoped_roster_policy="on")
    for f in d.families:
        assert not SCOPED & set(_names(f))


def test_an_unknown_policy_is_refused():
    with pytest.raises(ValueError):
        route(BMAYDIS, [], lineage_taxids_resolver=_lineage, scoped_roster_policy="yes")


def test_a_family_without_scoped_genes_is_the_same_object():
    f = Family(key=FamilyKey("X", "Y"), vocabulary_type="enum", idiomorph_values=["A"],
               idiomorph_pattern=None, genes=[{"name": "g", "role": "core_MAT"}],
               taxonomic_scope=[1])
    assert scope_roster_genes(f, set()) is f


def _query_names(path):
    return [line[1:] for line in path.read_text().splitlines() if line.startswith(">")]


def test_the_query_set_holds_the_scoped_proteins_only_when_the_routed_roster_does(tmp_path):
    db = Path("db")
    _, inside = _routed_mat(BMAYDIS, "on")
    _, outside = _routed_mat(AUREOBASIDIUM, "on")
    keys = {inside.key}
    default = _query_names(build_reference_fasta(db, tmp_path / "a.faa", family_keys=keys))
    out = _query_names(build_reference_fasta(db, tmp_path / "b.faa", family_keys=keys,
                                             families=[outside]))
    ins = _query_names(build_reference_fasta(db, tmp_path / "c.faa", family_keys=keys,
                                             families=[inside]))
    assert default == out, "outside the scope the query set is the unscoped one"
    assert not any(n.split("|")[2] in SCOPED for n in default)
    added = sorted(set(ins) - set(default))
    assert added == ["5016_c4_MAT_MAT1-2|gene1|GAP1", "5016_c4_MAT_MAT1-2|gene2|BGL1",
                     "5016_c5_MAT_MAT1-1|gene1|GAP1", "5016_c5_MAT_MAT1-1|gene2|BGL1"]
    searchable = searchable_genes_by_family(tmp_path / "c.faa", load_record_families(db))
    assert SCOPED <= searchable[inside.key]


def test_the_records_carry_the_deposit_proteins():
    """GAP1/BGL1 of the C4 (AF027687.1) and C5 (AF029913.1) deposits, 714 and 870 aa."""
    for rec in ("5016_c4_MAT_MAT1-2", "5016_c5_MAT_MAT1-1"):
        faa = (Path("db/Ascomycota/Pleosporales") / rec / "proteins.faa").read_text()
        seqs = {c.split("\n", 1)[0].split("|")[2].split("=")[1]: c.split("\n", 1)[1].replace("\n", "")
                for c in faa.split(">")[1:]}
        assert len(seqs["GAP1"]) == 714 and seqs["GAP1"].startswith("MQQEPAETMTLEQ")
        assert len(seqs["BGL1"]) == 870 and seqs["BGL1"].startswith("MASLLFLSLSGLAAG")


def _hit(gene, role, start):
    return SearchHit(family_key=FamilyKey("Ascomycota", "MAT"), gene_name=gene, role=role,
                     contig="c1", start=start, end=start + 900, strand="+", identity=90.0,
                     reference_record_id="r", method="tblastn_genome")


def _fraction(family, genes):
    roles = {"MAT1-1-1": "core_MAT"}
    cluster = GeneCluster("c1", 1, 40_000,
                          [_hit(g, roles.get(g, "flanking_conserved"), 1 + 5000 * i)
                           for i, g in enumerate(genes)])
    return score_cluster(cluster, [family])[0]


def test_a_found_slot_flank_raises_the_score_and_an_absent_one_costs_nothing():
    _, unscoped = _routed_mat(AUREOBASIDIUM, "on")
    _, scoped = _routed_mat(BMAYDIS, "on")
    flanked = ["MAT1-1-1", "APN2", "COX13"]
    assert _fraction(scoped, flanked).fraction_found == _fraction(unscoped, flanked).fraction_found == 0.75
    lone = _fraction(scoped, ["MAT1-1-1", "GAP1", "BGL1"])
    assert lone.fraction_found == 0.5 and {"GAP1", "BGL1"} <= set(lone.genes_found)
    assert _fraction(scoped, ["MAT1-1-1", "APN2", "COX13", "GAP1"]).fraction_found == 0.8


def test_a_gene_min_identity_drops_only_that_genes_weak_hits():
    """BGL1 is a beta-glucosidase family: a 27.7% fragment made spurious high calls."""
    from MATPredict.detect.search import drop_below_gene_min_identity
    _, scoped = _routed_mat(BMAYDIS, "on")
    _, unscoped = _routed_mat(AUREOBASIDIUM, "on")
    assert {g["name"]: g.get("min_identity") for g in scoped.genes if g["name"] in SCOPED} == {"GAP1": 50, "BGL1": 50}
    weak_bgl1 = _hit("BGL1", "flanking_conserved", 1)
    weak_bgl1 = type(weak_bgl1)(**{**weak_bgl1.__dict__, "identity": 27.7})
    strong_bgl1 = _hit("BGL1", "flanking_conserved", 5000)
    weak_core = type(weak_bgl1)(**{**_hit("MAT1-1-1", "core_MAT", 9000).__dict__, "identity": 20.0})
    kept = drop_below_gene_min_identity([weak_bgl1, strong_bgl1, weak_core], [scoped])
    assert kept == [strong_bgl1, weak_core]
    hits = [weak_core]
    assert drop_below_gene_min_identity(hits, [unscoped]) is hits, "no floor declared: untouched"


def test_a_failed_scope_lookup_is_recorded_not_silent():
    """A transient lookup error must not turn a policy-on run into an unmarked main run."""
    calls = {"n": 0}

    def flaky(taxid):
        calls["n"] += 1
        if calls["n"] > 1:  # routing's own lookup succeeds, the scope lookup fails
            raise TimeoutError("E-utilities")
        return _lineage(taxid)

    d = route(BMAYDIS, load_all_families(Path("db")), lineage_taxids_resolver=flaky,
              scoped_roster_policy="on")
    mat = next(f for f in d.families if f.key == FamilyKey("Ascomycota", "MAT"))
    assert not SCOPED & set(_names(mat))
    assert d.scoped_roster.startswith("SKIPPED") and "TimeoutError" in d.scoped_roster
    assert "scoped roster skipped" in d.routing_error


def test_the_decision_says_what_scoping_did():
    assert _routed_mat(BMAYDIS, "on")[0].scoped_roster == "on: BGL1,GAP1 searched (genome inside scope)"
    assert _routed_mat(AUREOBASIDIUM, "on")[0].scoped_roster == "on: BGL1,GAP1 not searched (genome outside scope)"
    assert _routed_mat(BMAYDIS, "off")[0].scoped_roster == "off: BGL1,GAP1 not searched"
    assert _routed_mat(BMAYDIS, "off")[0].routing_error is None


def test_records_carrying_a_scoped_gene_lie_inside_its_scope():
    from MATPredict.detect.scope_audit import audit_gene_scopes
    fams = load_all_families(Path("db"))
    assert audit_gene_scopes(Path("db"), fams, lineage_taxids_resolver=_lineage) == []
    outside = audit_gene_scopes(Path("db"), fams, lineage_taxids_resolver=lambda t: [4890, DOTHIDEALES])
    assert sorted((r, g) for r, g, _, _ in outside) == [
        ("5016_c4_MAT_MAT1-2", "BGL1"), ("5016_c4_MAT_MAT1-2", "GAP1"),
        ("5016_c5_MAT_MAT1-1", "BGL1"), ("5016_c5_MAT_MAT1-1", "GAP1")]
