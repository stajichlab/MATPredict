"""Collapsing redundant reference proteins onto one canonical gene name.

Curator's ruling, 2026-09-21: "make it nonredundant and collapse MFa1,MFa2,MFa3
into MFa and same for MFalpha1,MFalpha2 into MFalpha - for future collapse by
hash and pick a common name."

WHY, measured. A 334-genome Tremellales run reported 3,801 loci, of which 3,210
were a single genomic interval wearing several gene names -- 927 of them one
~95 bp ORF matching MFa1/MFa2/MFa3 or MFalpha1/MFalpha2. The cause is not a
ranking problem: those proteins are BYTE-IDENTICAL in the curated database
(MFa1/MFa2/MFa3 = 42 aa, one sequence; MFalpha1/MFalpha2 = 38 aa, one sequence),
so identity, coverage, e-value and bitscore are all exactly tied and no
best-hit rule can choose between them. The redundancy has to go at the source.

MFalpha3 differs by one residue and, per the curator's ruling, is folded into
the same group anyway: the GROUP is the unit of biological interest, and a
consensus or per-group HMM is the intended future search strategy. Collapsing
is a change of NAME, not of query set -- MFalpha3's distinct sequence is still
emitted and still searched, it simply reports as MFalpha.
"""
from pathlib import Path

import pytest

from MATPredict.detect.family_registry import load_all_families
from MATPredict.detect.reference_fasta import (
    build_reference_fasta,
    redundant_gene_name_groups,
    searchable_genes_by_family,
)


def _fam(phylum, locus):
    return next(f for f in load_all_families(Path("db"))
                if f.key.phylum == phylum and f.key.locus_name == locus)


def test_the_roster_declares_the_collapsed_names():
    genes = {g["name"] for g in _fam("Basidiomycota", "MAT").genes}
    assert "MFa" in genes and "MFalpha" in genes
    assert not {"MFa1", "MFa2", "MFa3", "MFalpha1", "MFalpha2", "MFalpha3"} & genes
    assert genes == {"SXI1", "SXI2", "MFalpha", "MFa", "STE3"}


def test_aliases_map_back_to_the_canonical_name():
    fam = _fam("Basidiomycota", "MAT")
    assert fam.gene_aliases["MFa2"] == "MFa"
    assert fam.gene_aliases["MFalpha1"] == "MFalpha"
    assert fam.gene_aliases["MFalpha3"] == "MFalpha"
    # a gene with no aliases maps to itself, so callers need no special case
    assert fam.gene_aliases.get("STE3", "STE3") == "STE3"


def test_the_reference_fasta_emits_one_protein_per_distinct_sequence(tmp_path):
    """Dedup is on (record, canonical name, SEQUENCE) -- not one per group.

    MFa1/MFa2/MFa3 are one sequence, so MFa goes out once. MFalpha covers two
    distinct sequences (MFalpha1==MFalpha2, and MFalpha3 one residue apart), so
    MFalpha goes out twice. Collapsing must not cost sensitivity: both variants
    are still searched, they just report under one name.
    """
    out = build_reference_fasta(Path("db"), tmp_path / "ref.faa")
    names = [line[1:].split("|")[2] for line in out.read_text().splitlines()
             if line.startswith(">")]
    assert names.count("MFa") == 1, "MFa1/MFa2/MFa3 are one sequence: emit it once"
    assert names.count("MFalpha") == 2, "two distinct sequences, both still searchable"
    for gone in ("MFa1", "MFa2", "MFa3", "MFalpha1", "MFalpha2", "MFalpha3"):
        assert gone not in names


def test_searchable_genes_sees_the_canonical_names(tmp_path):
    from MATPredict.detect.family_registry import load_record_families
    out = build_reference_fasta(Path("db"), tmp_path / "ref.faa")
    searchable = searchable_genes_by_family(out, load_record_families(Path("db")))
    genes = searchable[_fam("Basidiomycota", "MAT").key]
    assert "MFa" in genes and "MFalpha" in genes
    assert not {"MFa1", "MFa2", "MFalpha1"} & genes


def test_the_alpha_roster_shrinks_to_what_is_really_distinguishable(tmp_path):
    """The point of the whole change: one ORF can no longer count as 3 genes.

    Before, finding the single alpha pheromone ORF produced
    genes_found=[MFalpha1,MFalpha2,MFalpha3] -> 3/5 = 0.6, over the 0.5 floor.
    """
    from MATPredict.detect.family_registry import expected_genes_for_idiomorph
    fam = _fam("Basidiomycota", "MAT")
    alpha = [g["name"] for g in expected_genes_for_idiomorph(fam, ["SXI1"])]
    assert "MFalpha" in alpha
    assert not {"MFalpha1", "MFalpha2", "MFalpha3"} & set(alpha)
    assert alpha.count("MFalpha") == 1, "one roster slot, however many sequences back it"


def test_no_identical_sequence_group_spans_two_canonical_names():
    """The general rule, not just the two known cases: if two curated proteins
    are byte-identical they must resolve to ONE gene name, or a future
    curation will silently re-create this bug."""
    offenders = redundant_gene_name_groups(Path("db"))
    assert offenders == [], (
        "identical curated sequences under different canonical gene names; "
        "give them one name with `aliases:` in order.yml: " + repr(offenders)
    )
