"""Polished-evidence vote: a polished idiomorph beats an unpolished one (#275).

Curator's decision, 2026-11-06 (issue #275), behind `--polished-vote-policy {off,on}` (default `off`, in no preset until the release-prep
pull request). The idiomorph vote (`idiomorph.idiomorph_candidates`) ranks each idiomorph by the best bitscore of its informative hits, so
a short unpolished fragment of one idiomorph's core gene can outscore a polished model of the other idiomorph's: in the assembly-rescue
pilot an unpolished 47-bit MAT1-2-1 window beat a polished 39.7-bit MAT1-1-3 model in *Trichoderma viride* (truth MAT1-1), and the call
was MAT1-2.

The rule: when the called idiomorph rests ONLY on unpolished evidence (none of its informative core genes has a polished model) and the
other idiomorph has at least one polished informative core gene, the polished idiomorph is called. Nothing else changes except what
follows from a changed call the way it does for `--partial-core-policy` (the report's `idiomorph_vote_note` records the decision, the tier
is at most `medium`, and the vote's margin is not reported because it describes the idiomorph that was set aside).

Scope: `Ascomycota:MAT` only, the family it was measured on. Measured on the v0.8.1 run: of 3,742 `Ascomycota:MAT` loci where both
idiomorphs have a live core gene, 8 (all medium) match.

Everything here is pure; the pipeline supplies the per-gene facts.
"""
from __future__ import annotations

from dataclasses import dataclass

from MATPredict.detect.family_registry import Family
from MATPredict.detect.idiomorph import _informative

POLICY_OFF = "off"
POLICY_ON = "on"
POLISHED_VOTE_POLICIES = (POLICY_OFF, POLICY_ON)
DEFAULT_POLISHED_VOTE_POLICY = POLICY_OFF

#: The one family the rule applies to (phylum, locus_name).
SCOPE_FAMILY = ("Ascomycota", "MAT")
#: What the report's `idiomorph_vote_note` carries for a call this rule changed (the same vocabulary as `partial_core`).
RULE = "polished_vote"
DECISION_FLIPPED = "flipped"


def check_policy(policy: str) -> str:
    if policy not in POLISHED_VOTE_POLICIES:
        raise ValueError(f"unknown polished-vote policy {policy!r}; expected one of {POLISHED_VOTE_POLICIES}")
    return policy


@dataclass(frozen=True)
class GeneStat:
    """One live gene of the locus (not `lost_to`): its polish status and its bitscore."""

    gene: str
    status: str | None
    bitscore: float | None


@dataclass(frozen=True)
class Decision:
    idiomorph: str
    note: dict


def _polished(stat: GeneStat) -> bool:
    return (stat.status or "").startswith("polished")


def decide(family: Family, family_key: tuple[str, str], called: str, stats: list[GeneStat], *, policy: str = POLICY_ON) -> Decision | None:
    """None when the call stands; otherwise the polished idiomorph that replaces it."""
    if policy == POLICY_OFF or family_key != SCOPE_FAMILY or family.vocabulary_type != "enum" or called == "undetermined":
        return None
    by_name = {s.gene: s for s in stats}
    voters: dict[str, list[GeneStat]] = {}
    for gene in family.genes:
        idiomorphs = gene.get("present_in_idiomorphs") or []
        if len(idiomorphs) != 1 or not _informative(gene) or gene["name"] not in by_name:
            continue
        voters.setdefault(idiomorphs[0], []).append(by_name[gene["name"]])
    mine = voters.get(called, [])
    if not mine or any(_polished(s) for s in mine):
        return None
    others = {name: [s for s in genes if _polished(s)] for name, genes in voters.items() if name != called}
    others = {name: genes for name, genes in others.items() if genes}
    if not others:
        return None
    best_name = max(others, key=lambda n: (max(s.bitscore or 0.0 for s in others[n]), n))
    best = max(others[best_name], key=lambda s: s.bitscore or 0.0)
    called_best = max(mine, key=lambda s: s.bitscore or 0.0)
    return Decision(best_name, {
        "rule": RULE, "decision": DECISION_FLIPPED, "idiomorph": best_name, "from": called,
        "set_aside_gene": called_best.gene, "set_aside_status": called_best.status, "set_aside_bitscore": called_best.bitscore,
        "polished_gene": best.gene, "polished_status": best.status, "polished_bitscore": best.bitscore,
    })
