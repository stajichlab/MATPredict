"""btbA marks the Mucoromycota MAT locus but must not vote on the idiomorph.

Curator ruling, J. Stajich, 2026-09-26: mark btbA `idiomorph_informative:
false`.

WHY, measured on the 294-genome early-diverging Mucoromycota scan: 33 calls
carry sexM (about 96% identity, 99% coverage) and no sexP, yet were labelled
Plus. All 33 are Rhizopus, and btbA is present in all 33. btbA was marked
Plus-only because the one curated record that carries it (CBS 346-36) is a
Plus strain. It is a long protein, so its bitscore (e.g. 871 in
GCA_000696915.1) outvoted the sexM bitscore (378). A flank gene outvoted the
core MAT gene. See docs/notes/2026-09-26_early-diverging-scan-and-chytrid-control.md.
"""
from pathlib import Path

from MATPredict.detect.family_registry import load_all_families
from MATPredict.detect.idiomorph import (
    assign_idiomorph, evidenced_idiomorphs, idiomorph_candidates,
)
from MATPredict.detect.search import SearchHit

DB = Path(__file__).resolve().parents[2] / "db"


def _family():
    return next(f for f in load_all_families(DB)
                if (f.key.phylum, f.key.locus_name) == ("Mucoromycota", "MAT"))


def _h(fam, gene, role, bits, identity):
    return SearchHit(fam.key, gene, role, "c1", 1000, 1300, "+",
                     identity, "r", "tblastn_genome", bitscore=bits)


def test_btbA_is_non_informative_in_the_curated_roster():
    btbA = next(g for g in _family().genes if g["name"] == "btbA")
    assert btbA.get("idiomorph_informative") is False
    assert btbA.get("optional") is True          # still an optional locus marker


def test_a_sexM_only_rhizopus_locus_is_called_minus_despite_btbA():
    """GCA_000696915.1: btbA bitscore 871, sexM 378, sexP a 32% cross-hit."""
    fam = _family()
    hits = [_h(fam, "btbA", "flanking_variable", 871, 98.3),
            _h(fam, "sexM", "core_MAT", 378, 95.7),
            _h(fam, "sexP", "core_MAT", 40, 32.1)]
    names = [h.gene_name for h in hits]
    assert assign_idiomorph(fam, names, hits) == "Minus"
    ranked = idiomorph_candidates(fam, names, hits)
    assert ranked[0] == {"idiomorph": "Minus", "score": 378.0}


def test_btbA_alone_evidences_no_idiomorph():
    fam = _family()
    assert evidenced_idiomorphs(fam, [_h(fam, "btbA", "flanking_variable", 871, 98.3)]) == set()
