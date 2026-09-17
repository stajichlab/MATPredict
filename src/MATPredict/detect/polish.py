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


STATUS_AGREE = "polished_agree"
STATUS_DISAGREE = "polished_disagree"
STATUS_SINGLE = "polished_single"
STATUS_UNPOLISHED = "unpolished"


@dataclass(frozen=True)
class PolishOutcome:
    status: str
    canonical: PolishModel | None
    exonerate_model: PolishModel | None
    miniprot_model: PolishModel | None


def classify(
    exonerate_model: PolishModel | None,
    miniprot_model: PolishModel | None,
    tolerance_bp: int = 10,
) -> PolishOutcome:
    """Classify a gene's two-tool polish outcome. exonerate_model is
    preferred as canonical when both models exist (more directly
    integrated, matches sub-project 1's curation conventions) --
    miniprot's model is retained for comparison either way, never
    discarded. Both None with no fallback is an error at the caller
    (an unpolished status needs the raw tblastn hit, which this
    function does not have -- the caller in pipeline.py builds
    STATUS_UNPOLISHED's canonical from the tblastn SearchHit directly,
    not through this function)."""
    if exonerate_model is not None and miniprot_model is not None:
        status = STATUS_AGREE if boundaries_agree(exonerate_model, miniprot_model, tolerance_bp) else STATUS_DISAGREE
        return PolishOutcome(status, exonerate_model, exonerate_model, miniprot_model)
    if exonerate_model is not None or miniprot_model is not None:
        canonical = exonerate_model or miniprot_model
        return PolishOutcome(STATUS_SINGLE, canonical, exonerate_model, miniprot_model)
    return PolishOutcome(STATUS_UNPOLISHED, None, None, None)
