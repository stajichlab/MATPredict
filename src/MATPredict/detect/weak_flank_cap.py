"""A call whose only flank support is weak is not `high` (`--weak-flank-cap-policy`, #129).

Curator's decisions: 2026-10-22 (report only for the clean run) and 2026-11-07 (opt-in rule, each affected locus reviewed), issue #129.

A `high` Ascomycota:MAT call can rest on a "flank gene" that is noise: an unpolished APN2 HSP at 5-6% coverage and about 30 bits
(*Exserohilum*), or an unsuperseded opposite-idiomorph cross-hit. The extra gene lifts a 0.5 tie to `high` through the flanking branch of
`tiering.assign_tier`. Measured on the v0.8.1 Ascomycota run: of 12,369 `high` `Ascomycota:MAT` loci, 10,655 have a supported flank gene,
1,689 have no flank gene at all, and 25 have only weak ones.

With the policy `on`, a `high` `Ascomycota:MAT` call (idiomorph MAT1-1 or MAT1-2) that carries at least one live flank gene and in which
EVERY live flank gene is weak is lowered to `medium`. A flank gene is supported when its model is polished, or when it reaches both
`MIN_IDENTITY` identity and `MIN_COVERAGE` coverage; weak is neither. Only the tier changes; the idiomorph, the locus, the coordinates and
the class never do. The report's `weak_flank_cap` names the weak genes. A call with no flank gene is left alone (nothing weak lifted it).
"""
from __future__ import annotations

from dataclasses import dataclass

POLICY_OFF = "off"
POLICY_ON = "on"
WEAK_FLANK_CAP_POLICIES = (POLICY_OFF, POLICY_ON)
DEFAULT_WEAK_FLANK_CAP_POLICY = POLICY_OFF

SCOPE_FAMILY = ("Ascomycota", "MAT")
CALLED_IDIOMORPHS = ("MAT1-1", "MAT1-2")
FLANK_ROLE_PREFIX = "flanking"
#: A flank gene is supported by a polished model, or by this identity AND this coverage (percent).
MIN_IDENTITY = 50.0
MIN_COVERAGE = 30.0


def check_policy(policy: str) -> str:
    if policy not in WEAK_FLANK_CAP_POLICIES:
        raise ValueError(f"unknown weak-flank-cap policy {policy!r}; expected one of {WEAK_FLANK_CAP_POLICIES}")
    return policy


@dataclass(frozen=True)
class FlankStat:
    """What the rule needs to know of one live gene (not `lost_to`)."""

    gene: str
    role: str | None
    status: str | None
    identity: float | None
    coverage: float | None


def _supported(s: FlankStat) -> bool:
    return str(s.status or "").startswith("polished") or ((s.identity or 0.0) >= MIN_IDENTITY and (s.coverage or 0.0) >= MIN_COVERAGE)


def evaluate(confidence: str, idiomorph: str | None, family: tuple[str, str], stats: list[FlankStat], *, policy: str = POLICY_ON) -> dict | None:
    """None when the call keeps its tier; otherwise the note saying which flank genes were weak."""
    if policy == POLICY_OFF or confidence != "high":
        return None
    if family != SCOPE_FAMILY or idiomorph not in CALLED_IDIOMORPHS:
        return None
    flanks = [s for s in stats if str(s.role or "").startswith(FLANK_ROLE_PREFIX)]
    if not flanks or any(_supported(s) for s in flanks):
        return None
    return {
        "rule": "weak_flank_cap", "from": "high", "to": "medium",
        "weak_flank_genes": [{"gene": s.gene, "status": s.status, "identity": s.identity, "coverage": s.coverage} for s in flanks],
        "supported_if": f"polished, or identity >= {MIN_IDENTITY:g} and coverage >= {MIN_COVERAGE:g}",
    }
