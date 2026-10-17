"""A lone, strong core MAT gene as a low-confidence locus (`--strong-single-gene-policy`).

Curator's decisions 2026-10-17 (docs/notes/2026-10-17_single-gene-rule.md).
In a family that is otherwise unreported in a genome, ONE modelled core gene at
high identity to its curated reference, with no flank and no second gene in its
cluster, is reported as class `idiomorph_gene_only` at confidence `low`, flagged
`strong_single_gene`, instead of being dropped by the three gates in series
(evidence floor `min_hits=2`, the modelled-gene bar, the fraction floor).
Measured case: three Parastagonospora genomes (GCA_003497985.1, GCA_003503065.1,
GCA_022559985.1) carry the true gene at 99.4, 99.4 and 85.6% while HMG paralogs
in the same genomes sit at 40-56%.

OFF by default. Independent of `--idiomorph-policy` and `--partial-core-policy`.
It never changes a family that already has a report: it runs after the strict,
relaxed and split-locus passes and looks only at families with no call.

Parameters and why:

* identity gate 85 (curator). The identity is the MODELLED core gene's
  `GeneEvidence.identity` (the polished model's identity to its reference
  protein), the same quantity `split_locus` gates on; the raw tblastn hit
  identity of the cluster is recorded beside it as `cluster_identity`. A model
  is the evidence that stands for a locus, and a hit is only a localisation.
* coverage 90: the model must cover >= 90% of the length of the reference
  protein it was built from (`GeneEvidence.coverage`, a length ratio; None,
  i.e. unknown reference length, fails). "Near-full-length" is what separates
  a real gene from a conserved-domain fragment of a paralog; it is the same
  full-length bar `partial_core` uses (>= 90%).
* role `core_MAT`, status modelled (polished), idiomorph determined: a
  `flanking_*` gene, an unmodelled hit, or an `undetermined` idiomorph never
  qualifies.
* lone: the cluster's evidence holds no gene other than the core gene, except
  a hit that lost an idiomorph cross-match to it (`lost_to`).

When two qualifying loci of OPPOSITE idiomorph exist in one family, both are
reported (the best per idiomorph); the genome's `two_idiomorphs` statement is
then written for that family even if the roster does not enable it, because the
curator ruled such a genome is flagged, not resolved, and not called homothallic.
"""
from __future__ import annotations

from dataclasses import dataclass

from MATPredict.detect.split_locus import _is_modelled

POLICY_OFF = "off"
POLICY_ON = "on"
STRONG_SINGLE_GENE_POLICIES = (POLICY_OFF, POLICY_ON)
DEFAULT_STRONG_SINGLE_GENE_POLICY = POLICY_OFF

_UNDETERMINED = {None, "", "undetermined"}


@dataclass(frozen=True)
class StrongSingleGeneParams:
    min_core_identity: float = 85.0
    min_coverage: float = 90.0


DEFAULT_STRONG_SINGLE_GENE = StrongSingleGeneParams()


def evaluate_strong_single_gene(
    candidate, *, params: StrongSingleGeneParams = DEFAULT_STRONG_SINGLE_GENE,
) -> dict | None:
    """The `strong_single_gene` report block when `candidate` qualifies, else None.

    `candidate` is a built DetectionResult for one cluster."""
    if candidate.idiomorph in _UNDETERMINED:
        return None
    live = [e for e in candidate.gene_evidence if e.lost_to is None]
    if len(live) != 1:
        return None
    core = live[0]
    if core.role != "core_MAT" or not _is_modelled(core):
        return None
    if core.identity is None or core.identity < params.min_core_identity:
        return None
    if core.coverage is None or core.coverage < params.min_coverage:
        return None
    return {
        "core_gene": core.gene_name,
        "core_contig": core.contig,
        "core_start": core.start,
        "core_end": core.end,
        "core_identity": round(core.identity, 3),
        "core_coverage": core.coverage,
        "reference_record": core.reference_record_id,
        "min_core_identity": params.min_core_identity,
        "min_coverage": params.min_coverage,
    }
