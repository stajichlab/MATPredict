"""Cross-tool gene-model comparison for the localize-then-polish search
revision -- see docs/superpowers/specs/2026-09-17-mat-detection-search-localization-design.md."""
from __future__ import annotations

from dataclasses import dataclass

from MATPredict.detect.family_registry import FamilyKey


@dataclass(frozen=True)
class ExonSpan:
    start: int
    end: int


@dataclass(frozen=True)
class PolishModel:
    gene_name: str
    family_key: FamilyKey
    role: str
    contig: str
    start: int
    end: int
    strand: str
    exons: list[ExonSpan]
    identity: float
    reference_record_id: str
    method: str


def boundaries_agree(a: PolishModel, b: PolishModel, tolerance_bp: int = 10) -> bool:
    """Same exon count, and each corresponding exon's start/end within tolerance_bp."""
    if len(a.exons) != len(b.exons):
        return False
    return all(
        abs(ea.start - eb.start) <= tolerance_bp and abs(ea.end - eb.end) <= tolerance_bp
        for ea, eb in zip(a.exons, b.exons)
    )
