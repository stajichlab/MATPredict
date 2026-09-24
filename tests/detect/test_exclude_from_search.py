"""`exclude_from_search: true` -- a gene the records keep but no run queries.

Curator's ruling, 2026-09-21: "don't necessarily delete RPO41 from the library
but keep it off the search list for now."

Distinct from the two flags that already exist:

* `optional: true` -- searched and reported, but out of `fraction_found` and
  out of the core requirement. The biology does not require it.
* `genes_not_searchable` -- a run HELD no reference protein for it, so it could
  not have been found. That is a fact about the reference set, not a choice.

This is a deliberate choice to stop querying a protein whose curation is sound.
Measured reason on Cryptococcus JEC21: STE20 (a p21-activated kinase) produced
187 of the 272 genome-wide tblastn HSPs, on all 14 contigs, at 31.7% median
identity, and every one of those fed the polish stage -- a run went from 8 s to
535 s. RPO41 was excluded for a different reason: it is well behaved (14 HSPs,
one contig, 98.7% median) but is not established as universal outside
Cryptococcus.

The record keeps the gene, so the curation is not lost and the gene can be
restored by deleting one flag. The roster keeps declaring it, so
`validate_gene_vocabulary` still recognises the curated name.
"""
from pathlib import Path

from MATPredict.detect.family_registry import load_all_families, load_record_families
from MATPredict.detect.reference_fasta import (
    build_reference_fasta,
    searchable_genes_by_family,
)

EXCLUDED = {"STE20", "RPO41", "UAP1", "NOG2"}


def _mat():
    return next(f for f in load_all_families(Path("db"))
                if f.key.phylum == "Basidiomycota" and f.key.locus_name == "MAT")


def test_the_roster_still_declares_the_excluded_genes():
    """They must stay declared, or the schema guard would call the curated
    records' own gene names orphaned."""
    genes = {g["name"] for g in _mat().genes}
    assert EXCLUDED <= genes


def test_the_records_still_hold_their_proteins():
    for rec in ("40410_jec21_MAT_alpha", "40410_jec20_MAT_a"):
        faa = Path("db/Basidiomycota/Tremellales") / rec / "proteins.faa"
        names = {line.split("|")[2].split("=")[1]
                 for line in faa.read_text().splitlines() if line.startswith(">")}
        assert EXCLUDED <= names, f"{rec} must keep the curated proteins"


def test_the_reference_fasta_omits_them(tmp_path):
    out = build_reference_fasta(Path("db"), tmp_path / "ref.faa")
    names = {line[1:].split("|")[2] for line in out.read_text().splitlines()
             if line.startswith(">")}
    assert not (EXCLUDED & names), "excluded genes must not be queried"
    # the ones that stay
    assert {"FAO1", "PAN6", "MFalpha", "MFa", "STE3", "SXI1", "RPL39"} <= names


def test_they_are_reported_as_not_searchable(tmp_path):
    """Falling out of the reference FASTA puts them in `genes_not_searchable`
    by the existing route, so they leave the denominator and the core
    requirement without any new special case."""
    out = build_reference_fasta(Path("db"), tmp_path / "ref.faa")
    searchable = searchable_genes_by_family(out, load_record_families(Path("db")))
    assert not (EXCLUDED & searchable[_mat().key])
