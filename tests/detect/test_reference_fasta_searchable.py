# tests/detect/test_reference_fasta_searchable.py
"""Which genes a run could have found, read from the FASTA it searched with.

Derived from the written reference FASTA rather than from the database, so the
query set and the scoring denominator cannot disagree -- the same guarantee
`build_reference_fasta` already makes for hit attribution. Phylum routing
narrows this file (181 proteins to 19 for Mucoromycota), so a gene can have a
reference in db/ and still not be searchable in a given run.
"""
from __future__ import annotations

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.reference_fasta import searchable_genes_by_family

MAT = FamilyKey("Mucoromycota", "MAT")
ALPHA = FamilyKey("Basidiomycota", "Aalpha")


def _fasta(tmp_path, text):
    path = tmp_path / "_reference.faa"
    path.write_text(text)
    return path


def test_gene_names_are_grouped_by_the_family_their_record_belongs_to(tmp_path):
    fasta = _fasta(tmp_path, ">recA|gene0|tptA\nMKV\n>recA|gene1|sexP\nMKV\n>recB|gene0|Z\nMKV\n")
    result = searchable_genes_by_family(fasta, {"recA": MAT, "recB": ALPHA})
    assert result == {MAT: {"tptA", "sexP"}, ALPHA: {"Z"}}


def test_one_gene_contributed_by_several_records_appears_once(tmp_path):
    fasta = _fasta(tmp_path, ">recA|gene0|tptA\nMKV\n>recB|gene0|tptA\nMKV\n")
    result = searchable_genes_by_family(fasta, {"recA": MAT, "recB": MAT})
    assert result == {MAT: {"tptA"}}


def test_a_record_with_no_resolvable_family_is_skipped(tmp_path):
    # Such a protein could not be attributed to a family by search._attribute
    # either, so it cannot make any gene findable.
    fasta = _fasta(tmp_path, ">recA|gene0|tptA\nMKV\n>orphan|gene0|mystery\nMKV\n")
    result = searchable_genes_by_family(fasta, {"recA": MAT})
    assert result == {MAT: {"tptA"}}


def test_an_empty_reference_yields_no_searchable_genes(tmp_path):
    assert searchable_genes_by_family(_fasta(tmp_path, ""), {"recA": MAT}) == {}


def test_a_missing_reference_file_yields_no_information_and_warns(tmp_path, caplog):
    # A caller that stubs out the search never writes this file. Absence must
    # degrade to "no searchability information" -- which scoring already reads
    # as "keep the whole roster", the prior behaviour -- rather than aborting
    # the run. It warns, because in a real run the file is always written
    # before the search and its absence would mean something else went wrong.
    import logging

    with caplog.at_level(logging.WARNING):
        result = searchable_genes_by_family(tmp_path / "absent.faa", {"recA": MAT})
    assert result == {}
    assert "absent.faa" in caplog.text


def test_the_real_mucoromycota_reference_omits_algA_and_glrA(tmp_path):
    # The defect this exists to expose, asserted against the live database:
    # two of the seven genes in the Mucoromycota MAT roster have no reference
    # protein, so no genome can ever match them.
    from pathlib import Path

    from MATPredict.detect.family_registry import load_record_families
    from MATPredict.detect.reference_fasta import build_reference_fasta

    db_root = Path(__file__).parents[2] / "db"
    out = tmp_path / "ref.faa"
    build_reference_fasta(db_root, out, family_keys={MAT})
    searchable = searchable_genes_by_family(out, load_record_families(db_root))
    assert searchable[MAT] == {"tptA", "rnhA", "sexP", "sexM", "btbA"}
    assert "algA" not in searchable[MAT]
    assert "glrA" not in searchable[MAT]
