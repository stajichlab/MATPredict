"""A close idiomorph vote between short fragments is reported `undetermined` (`--short-fragment-vote-policy`, #138).

Curator's decision, 2026-11-07 (issue #138): "fragments only". Behind `--short-fragment-vote-policy {off,undetermined}` (default `off`).

The idiomorph of a locus is the winner of a bit-score vote between its candidate idiomorphs. When the winner beats the runner-up by less
than the family's minimum margin (`min_idiomorph_margin`, 5 bits) AND the best core gene of the locus covers less than `MAX_COVERAGE` of its
reference, the vote rests on short fragments a few bits apart and the label is a coin toss (*Aureobasidium* is homothallic, so both
labels are wrong; two *Xylaria* genomes swap labels on 30-50 aa HMG fragments). The call then reports `undetermined`; the vote's winner,
the margin and the best core gene move into `idiomorph_vote_note`. The locus, its coordinates and its tier are not changed here.

Measured on the v0.8.1 Ascomycota run: 158 `Ascomycota:MAT` loci have a multi-candidate margin under 5 bits; 108 have a best core gene
under 60% coverage (these are affected) and 50 have a fuller one (kept). A vote with a fuller-length gene is a real call, however close.
"""
from __future__ import annotations

POLICY_OFF = "off"
POLICY_UNDETERMINED = "undetermined"
SHORT_FRAGMENT_VOTE_POLICIES = (POLICY_OFF, POLICY_UNDETERMINED)
DEFAULT_SHORT_FRAGMENT_VOTE_POLICY = POLICY_OFF

#: Percent of the reference the best core gene must cover for a close vote to stand.
MAX_COVERAGE = 60.0


def check_policy(policy: str) -> str:
    if policy not in SHORT_FRAGMENT_VOTE_POLICIES:
        raise ValueError(f"unknown short-fragment-vote policy {policy!r}; expected one of {SHORT_FRAGMENT_VOTE_POLICIES}")
    return policy


def decide(policy: str, idiomorph: str | None, margin: float | None, n_candidates: int, min_margin: float,
           core_coverages: list[float | None]) -> dict | None:
    """The `idiomorph_vote_note` when the rule makes the call `undetermined`, else None."""
    if policy == POLICY_OFF or not idiomorph or idiomorph == "undetermined":
        return None
    if margin is None or n_candidates < 2 or margin >= min_margin:
        return None
    covs = [c for c in core_coverages if c is not None]
    if not covs or max(covs) >= MAX_COVERAGE:
        return None
    return {
        "rule": "short_fragment_vote", "decision": "undetermined", "vote_idiomorph": idiomorph, "margin": round(margin, 3),
        "min_margin": min_margin, "best_core_coverage": round(max(covs), 1), "max_coverage": MAX_COVERAGE,
    }
