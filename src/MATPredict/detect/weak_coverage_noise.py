"""An unpolished gene that covers little of its reference is noise, not a failed model (`--weak-coverage-noise-policy`, #149).

Curator's decision, 2026-11-07 (issue #149): "weak by coverage", opt-in. Behind `--weak-coverage-noise-policy {off,on}` (default `off`).

The pipeline treats a gene whose every hit scores under the family's noise floor (`flank_carried_min_bitscore`, 39 bits) as noise: it does
not cap the call at `medium`. A gene at or above the floor that no tool could model does cap it. So the same weak fragment gives a
different tier depending only on which side of 39 bits it falls, and a better reference that lifts a fragment from 36 to 42 bits can LOWER
a call's confidence although nothing about the evidence got worse (three Serinales A loci, #147). Measured on the v0.8.1 Ascomycota run: of
2,885 `medium` `Ascomycota:MAT` loci, 681 have their strongest counted unpolished gene at 39-45 bits.

With the policy `on`, a gene whose every hit covers under `MIN_COVERAGE` percent of its reference is also noise, whatever its bitscore. The
line is the same 30% as the weak-flank cap (#129). A hit with unknown coverage never qualifies. Only the tier can change (`medium` to
`high`, never the reverse); the idiomorph, the locus and the coordinates never do.
"""
from __future__ import annotations

POLICY_OFF = "off"
POLICY_ON = "on"
WEAK_COVERAGE_NOISE_POLICIES = (POLICY_OFF, POLICY_ON)
DEFAULT_WEAK_COVERAGE_NOISE_POLICY = POLICY_OFF

#: A gene whose every hit covers less than this percent of its reference is noise.
MIN_COVERAGE = 30.0


def check_policy(policy: str) -> str:
    if policy not in WEAK_COVERAGE_NOISE_POLICIES:
        raise ValueError(f"unknown weak-coverage-noise policy {policy!r}; expected one of {WEAK_COVERAGE_NOISE_POLICIES}")
    return policy


def low_coverage_genes(hits, family_key, min_coverage: float = MIN_COVERAGE) -> set[str]:
    """This family's genes whose EVERY live hit in `hits` has a known coverage under `min_coverage` percent."""
    best: dict[str, float | None] = {}
    for h in hits:
        if h.family_key != family_key or getattr(h, "superseded_by", None) is not None:
            continue
        cov = h.coverage
        if cov is None:
            best[h.gene_name] = None
        elif h.gene_name not in best:
            best[h.gene_name] = float(cov)
        elif best[h.gene_name] is not None:
            best[h.gene_name] = max(best[h.gene_name], float(cov))
    return {g for g, c in best.items() if c is not None and c < min_coverage}
