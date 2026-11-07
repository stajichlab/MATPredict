"""Report a single-candidate identity gap apart from the idiomorph vote margin (`--idiomorph-margin-policy`, #145).

Curator's decision, 2026-11-07 (issue #145), behind `--idiomorph-margin-policy {off,separate}` (default `off`).

`idiomorph_margin` normally is the bit-score gap by which the winning idiomorph beat the runner-up. For a locus with a SINGLE idiomorph
candidate there is no runner-up; the pipeline then falls back to the narrowest overlap resolution's identity gap (identity points between
two different genes on one HMG box). That number is a different quantity and is not a vote against the call. Measured on the v0.8.1
Ascomycota run: 6,684 of 18,107 `Ascomycota:MAT` loci carry such a fallback margin (4,984 high, 1,679 medium, 21 low), 986 of them under 5.
The tier cap for a thin margin already skips it (#137); this rule only fixes what the report says.

With the policy `separate` the fallback is reported as `idiomorph_identity_gap` and `idiomorph_margin` is `null`. Nothing else changes:
not the idiomorph, the locus, the coordinates, the tier or the class. A margin from a real vote (two or more candidates) or from the
profile-HMM classifier is never touched.
"""
from __future__ import annotations

from dataclasses import replace

POLICY_OFF = "off"
POLICY_SEPARATE = "separate"
IDIOMORPH_MARGIN_POLICIES = (POLICY_OFF, POLICY_SEPARATE)
DEFAULT_IDIOMORPH_MARGIN_POLICY = POLICY_OFF


def check_policy(policy: str) -> str:
    if policy not in IDIOMORPH_MARGIN_POLICIES:
        raise ValueError(f"unknown idiomorph-margin policy {policy!r}; expected one of {IDIOMORPH_MARGIN_POLICIES}")
    return policy


def is_identity_gap(r) -> bool:
    """The reported margin of call `r` is the single-candidate identity-gap fallback, not a vote margin or a classifier margin."""
    return (r.idiomorph_margin is not None and len(r.idiomorph_candidates or []) <= 1
            and not r.idiomorph_classifier and bool(r.idiomorph_resolutions))


def apply(results: list, policy: str) -> list:
    """`results` with each fallback margin moved to `idiomorph_identity_gap`; `off` returns them untouched."""
    if check_policy(policy) == POLICY_OFF:
        return results
    return [replace(r, idiomorph_margin=None, idiomorph_identity_gap=r.idiomorph_margin) if is_identity_gap(r) else r for r in results]
