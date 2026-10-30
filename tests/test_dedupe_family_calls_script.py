"""scripts/dedupe_family_calls.py: one locus called by two families keeps the better-scoring family."""
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("dedupe_family_calls", Path(__file__).resolve().parent.parent / "scripts" / "dedupe_family_calls.py")
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)


def row(fam, start="20", end="100", cls="idiomorph_gene_only", conf="high", polished="2", genes="A|B", genome="g1", contig="c1"):
    return dict(genome=genome, contig=contig, start=start, end=end, family_called=fam, locus_class=cls, confidence=conf, polished_genes=polished, genes_found=genes)


def test_better_class_keeps_the_locus():
    rows = mod.mark_duplicates([row("Ascomycota:MATsc", cls="mat_locus"), row("Ascomycota:MTL")])
    assert [r["duplicate_family_call"] for r in rows] == ["", "true"]
    assert {r["kept_family"] for r in rows} == {"Ascomycota:MATsc"}


def test_tie_goes_to_the_family_with_more_loci_in_the_genome():
    rows = mod.mark_duplicates([row("Ascomycota:MATsc"), row("Ascomycota:MTL"), row("Ascomycota:MTL", start="500", end="900")])
    dup = [r for r in rows if r["duplicate_family_call"]]
    assert len(dup) == 1 and dup[0]["family_called"] == "Ascomycota:MATsc"


def test_different_coordinates_or_one_family_are_not_duplicates():
    rows = mod.mark_duplicates([row("Ascomycota:MATsc"), row("Ascomycota:MTL", start="30"), row("Ascomycota:MATsc", start="900", end="999", genome="g2")])
    assert all(r["duplicate_family_call"] == "" for r in rows)


def test_overlap_mode_groups_nested_and_partly_overlapping_loci():
    # the same locus with different boundaries: one nested in the other (MATsc 100-1000 against MTL 200-900), and a short MTL locus that only just touches
    rows = mod.mark_duplicates([row("Ascomycota:MATsc", start="100", end="1000", cls="mat_locus"), row("Ascomycota:MTL", start="200", end="900"), row("Ascomycota:MTL", start="990", end="5000", cls="mat_locus")], match="overlap", min_overlap=0.9)
    assert [r["duplicate_family_call"] for r in rows] == ["", "true", ""]      # the touching locus overlaps 11 bp of 911: a second locus, not a duplicate


def test_overlap_mode_is_off_by_default_and_needs_different_families():
    nested = [row("Ascomycota:MATsc", start="100", end="1000"), row("Ascomycota:MTL", start="200", end="900")]
    assert all(r["duplicate_family_call"] == "" for r in mod.mark_duplicates([dict(r) for r in nested]))
    same_family = mod.mark_duplicates([row("Ascomycota:MATsc", start="100", end="1000"), row("Ascomycota:MATsc", start="200", end="900")], match="overlap")
    assert all(r["duplicate_family_call"] == "" for r in same_family)


def test_overlap_mode_requires_the_same_contig_and_genome():
    rows = mod.mark_duplicates([row("Ascomycota:MATsc", start="100", end="1000"), row("Ascomycota:MTL", start="200", end="900", contig="c2"), row("Ascomycota:MTL", start="200", end="900", genome="g2")], match="overlap")
    assert all(r["duplicate_family_call"] == "" for r in rows)
