"""Idiomorph assignment from which idiomorph-restricted genes were found."""
from __future__ import annotations

from MATPredict.detect.family_registry import Family


def assign_idiomorph(family: Family, genes_found: list[str]) -> str:
    if family.vocabulary_type != "enum":
        return "undetermined"  # pattern (multiallelic) families: allele number isn't callable by homology alone

    indicated: set[str] = set()
    for gene in family.genes:
        if gene["name"] in genes_found:
            indicated.update(gene.get("present_in_idiomorphs", []))

    if len(indicated) == 1:
        return next(iter(indicated))
    return "undetermined"
