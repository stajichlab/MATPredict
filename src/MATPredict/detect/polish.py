"""Cross-tool gene-model comparison for the localize-then-polish search
revision -- see docs/superpowers/specs/2026-09-17-mat-detection-search-localization-design.md.

Five report-facing `GeneEvidence.status` values are defined here:

* `STATUS_AGREE` -- both `exonerate --refine` and `miniprot` produced a
  model and the two agree (same exon count, boundaries within tolerance).
* `STATUS_DISAGREE` -- both tools produced a model but they disagree; the
  non-canonical tool's model is retained on `GeneEvidence.alternate_model`
  for a human reviewer to inspect.
* `STATUS_SINGLE` -- only one of the two tools produced a model.
* `STATUS_UNPOLISHED` -- the gene WAS a polish candidate (localized by
  `tblastn`, or a rescue target for a family's own missing `core_MAT`
  gene) and was sent through both polishing tools, but NEITHER tool could
  produce a usable model. This is a genuinely uncertain result: it caps
  the family's confidence tier at Medium (`_any_gene_unpolished` in
  `pipeline.py`, which reads the internal `PolishOutcome.status`, never
  `GeneEvidence.status`).
* `STATUS_NOT_POLISH_CANDIDATE` -- the gene was evidenced directly (e.g. a
  confident fast-path `diamond` hit, or an already-present core gene) and
  never entered the localize/polish pipeline at all -- no `PolishOutcome`
  was ever computed for it. This is a solid result with no implied
  uncertainty, and is reported so a curator does not confuse it with
  `STATUS_UNPOLISHED`'s genuinely-uncertain meaning.
"""
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
#: A gene evidenced directly (e.g. a confident fast-path diamond hit, or an
#: already-present core gene) that never entered the localize/polish
#: pipeline at all -- no PolishOutcome was ever computed for it. Distinct
#: from STATUS_UNPOLISHED, which is reserved for a gene that WAS localized
#: and sent through both polish tools but that neither tool could confirm
#: (see module docstring). This constant is purely report-facing: tiering's
#: `_any_gene_unpolished` in pipeline.py never reads it.
STATUS_NOT_POLISH_CANDIDATE = "not_polish_candidate"


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
