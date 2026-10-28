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
