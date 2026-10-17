"""Confidence cap for calls resting on partial or weak core evidence.

Curator's decision, 2026-10-17: cap `high` at `medium` (1) when the best core
evidence of the call is PARTIAL coverage and (2) when it is a WEAK-identity hit
with no flank gene. Behind `--confidence-cap-policy {off,partial_weak,weak_only}`
(default `off`); the criteria below were chosen from a replay on the v0.6.0
reports, argued in docs/notes/2026-10-17_confidence-cap.md.

What the data allowed (v0.6.0, 13,740 `high` loci):

* A bare "best core covers < 80% of its reference" would demote 49% of all
  `high` calls: a domain-only model against a long reference from another clade
  is ordinary (the same finding that shaped `partial_core`). On the BioSample
  truth set it demotes ~15-25 right calls per wrong one. It is therefore only
  applied where the partialness is demonstrably real: the gene lies within
  `EDGE_BP` of a contig end (an assembly truncation), and the call carries no
  full-length (>= `FULL_COVERAGE`) core gene of its own.
* The weak-identity criterion applies to calls with no flank gene, which in
  practice means `idiomorph_gene_only` loci (a `mat_locus` always has a flank).

The cap only ever lowers `high` to `medium`. It never raises, never touches
`medium`/`low`, never changes the locus, the idiomorph or the class, and is
idempotent with the narrow caps already in force (`partial_core`'s flipped /
undetermined calls, `cap_at_medium` on partial-strength calls): the later
caps are evaluated on whatever tier is left, so a call already at `medium`
is left alone and carries no note.

Everything here is pure; the pipeline supplies the per-gene facts.
"""
from __future__ import annotations

from dataclasses import dataclass

from MATPredict.detect.partial_core import EDGE_BP, FULL_COVERAGE, PARTIAL_COVERAGE

POLICY_OFF = "off"
POLICY_PARTIAL_WEAK = "partial_weak"
POLICY_WEAK_ONLY = "weak_only"
CONFIDENCE_CAP_POLICIES = (POLICY_OFF, POLICY_PARTIAL_WEAK, POLICY_WEAK_ONLY)
DEFAULT_CONFIDENCE_CAP_POLICY = POLICY_OFF

#: Percent identity below which the best core hit of a flank-less call is "weak".
#: The curator's "about 45". PROVISIONAL: the v0.6.0 identity distribution of
#: `high` gene-only loci is bimodal with its trough at 65-80, not at 45 (119 loci
#: below 45, 162 below 50, 245 below 65); see the note.
WEAK_IDENTITY = 45.0

CRITERION_PARTIAL = "partial_edge"
CRITERION_WEAK = "weak_identity"


def check_policy(policy: str) -> str:
    if policy not in CONFIDENCE_CAP_POLICIES:
        raise ValueError(f"unknown confidence-cap policy {policy!r}; expected one of {CONFIDENCE_CAP_POLICIES}")
    return policy


@dataclass(frozen=True)
class CapStat:
    """What the rule needs to know of one core gene of the CALLED allele (modelled or not, not `lost_to`)."""

    gene: str
    coverage: float | None  # percent of reference; None when unknown
    identity: float | None
    bitscore: float | None
    edge_distance: int | None  # bp to the nearer contig end; None when unknown


@dataclass(frozen=True)
class CapDecision:
    criteria: tuple[str, ...]
    note: dict


def evaluate(
    confidence: str,
    stats: list[CapStat],
    n_flank_genes: int,
    *,
    policy: str = POLICY_PARTIAL_WEAK,
    partial: float = PARTIAL_COVERAGE,
    full: float = FULL_COVERAGE,
    edge_bp: int = EDGE_BP,
    weak_identity: float = WEAK_IDENTITY,
) -> CapDecision | None:
    """None when the call keeps its confidence. Otherwise why it drops to `medium`.

    `stats` are the called allele's core genes; the best is the highest bitscore.
    `n_flank_genes` counts flanking genes (either flank role) found in the call.
    """
    if policy == POLICY_OFF or confidence != "high" or not stats:
        return None
    best = max(stats, key=lambda s: (s.bitscore or 0.0))
    criteria: list[str] = []
    if (policy == POLICY_PARTIAL_WEAK
            and best.coverage is not None and best.coverage < partial
            and not any(s.coverage is not None and s.coverage >= full for s in stats)
            and best.edge_distance is not None and best.edge_distance < edge_bp):
        criteria.append(CRITERION_PARTIAL)
    if n_flank_genes == 0 and best.identity is not None and best.identity < weak_identity:
        criteria.append(CRITERION_WEAK)
    if not criteria:
        return None
    return CapDecision(tuple(criteria), {
        "rule": policy, "from": "high", "to": "medium", "criteria": list(criteria),
        "best_core_gene": best.gene, "best_core_coverage": best.coverage, "best_core_identity": best.identity,
        "best_core_edge_distance": best.edge_distance, "flank_genes": n_flank_genes,
        "thresholds": {"partial_coverage": partial, "full_coverage": full, "edge_bp": edge_bp,
                       "weak_identity": weak_identity, "weak_identity_provisional": True},
    })
