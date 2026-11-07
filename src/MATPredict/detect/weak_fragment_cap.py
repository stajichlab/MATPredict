"""Weak-fragment tier cap: a call whose strongest core gene is a thin fragment is capped at `low` (#275).

Curator's decisions, 2026-11-06 (issue #275), behind `--weak-fragment-cap-policy {off,on}` (default `off`, in no preset until the
release-prep pull request):

* Scope: `Ascomycota:MAT` only, the family it was measured on. Other families keep their behaviour until someone measures them.
* Trigger: the STRONGEST core gene of the locus (a `MAT1-*` gene that is not `lost_to`, best by polish status -- both tools agree,
  then one tool or the tools disagree, then unpolished -- and then by bitscore) is not a both-tools-agree model and scores under
  `MAX_BITSCORE` bits.
* Effect: the tier is lowered to `low`; the idiomorph label, the locus and the class are never changed. A call already at `low` is
  untouched, and nothing is ever raised.

Measured on the 12,051 `Ascomycota:MAT` loci with a called idiomorph of the v0.8.1 run: 108 match, 40 of them already `low`, 68 would
change (23 medium, 45 high), 0.6% of the loci. A trigger of "no core gene polished" would have changed one locus; "identity under 40%"
would have changed 1,626. The reasoning and the alternatives are in docs/notes/2026-11-06_weak-fragment-cap.md.

Everything here is pure; the pipeline supplies the per-gene facts.
"""
from __future__ import annotations

from dataclasses import dataclass

from MATPredict.detect.polish import STATUS_AGREE

POLICY_OFF = "off"
POLICY_ON = "on"
WEAK_FRAGMENT_CAP_POLICIES = (POLICY_OFF, POLICY_ON)
DEFAULT_WEAK_FRAGMENT_CAP_POLICY = POLICY_OFF

#: The one family the rule applies to (phylum, locus_name), the idiomorph labels it applies to, and the prefix of its core genes.
SCOPE_FAMILY = ("Ascomycota", "MAT")
CALLED_IDIOMORPHS = ("MAT1-1", "MAT1-2")
CORE_GENE_PREFIX = "MAT1-"
#: Bitscore under which the strongest core gene is a "weak fragment" (the pipeline's own weak-fragment line: its noise floor is 39).
MAX_BITSCORE = 50.0


def check_policy(policy: str) -> str:
    if policy not in WEAK_FRAGMENT_CAP_POLICIES:
        raise ValueError(f"unknown weak-fragment-cap policy {policy!r}; expected one of {WEAK_FRAGMENT_CAP_POLICIES}")
    return policy


@dataclass(frozen=True)
class CoreStat:
    """What the rule needs to know of one live gene of the locus (not `lost_to`)."""

    gene: str
    status: str | None  # one of polish.STATUS_*
    bitscore: float | None


@dataclass(frozen=True)
class CapDecision:
    note: dict


def _strength(stat: CoreStat) -> tuple[int, float]:
    status = stat.status or ""
    rank = 2 if status == STATUS_AGREE else 1 if status.startswith("polished") else 0
    return (rank, stat.bitscore or 0.0)


def evaluate(
    confidence: str,
    idiomorph: str | None,
    family: tuple[str, str],
    stats: list[CoreStat],
    *,
    policy: str = POLICY_ON,
    max_bitscore: float | None = None,
) -> CapDecision | None:
    """None when the call keeps its confidence; otherwise why it drops to `low`."""
    if policy == POLICY_OFF or confidence not in ("medium", "high"):
        return None
    if family != SCOPE_FAMILY or idiomorph not in CALLED_IDIOMORPHS:
        return None
    limit = MAX_BITSCORE if max_bitscore is None else max_bitscore
    core = [s for s in stats if s.gene.startswith(CORE_GENE_PREFIX)]
    if not core:
        return None
    strongest = max(core, key=_strength)
    if strongest.status == STATUS_AGREE or (strongest.bitscore or 0.0) >= limit:
        return None
    return CapDecision({
        "rule": "weak_fragment_cap", "from": confidence, "to": "low",
        "strongest_core_gene": strongest.gene, "strongest_core_status": strongest.status,
        "strongest_core_bitscore": strongest.bitscore, "threshold_bitscore": limit,
    })
