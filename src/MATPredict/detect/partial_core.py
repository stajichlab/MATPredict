"""Partial-coverage vote rule: a truncated core gene must not beat a full-length
core gene of the other idiomorph.

Curator's decision, 2026-10-17. The idiomorph vote (`idiomorph.idiomorph_candidates`)
takes the best single bitscore per idiomorph. A 5'-truncated remnant of one
idiomorph's core gene (Botrytis cinerea: a 214-aa piece of MAT1-1-1, 61% of its
reference, 92.7% identical to its close reference) can then outscore the real,
full-length core gene of the other idiomorph (MAT1-2-1, 376 aa, 99.7% to the
curated record) because the reference sets are unequal and the bitscore does not
see how much of the gene is there. 17 BioProject PRJNA624742 genomes stated
MAT1-2 were called MAT1-1 (high) this way.

The rule, behind `--partial-core-policy coverage_profile` (default `off`):

1. TRIGGER (a coverage rule). The vote's winning idiomorph W has a best core hit
   whose model covers less than `PARTIAL_COVERAGE` of its reference protein, W
   has no full-length core gene, and the other idiomorph O has a full-length
   (>= `FULL_COVERAGE`) modelled core gene.
2. GUARD. If the partial hit lies within `EDGE_BP` of a contig end it is
   treated as an assembly artefact (a true gene cut by the contig); the old vote
   stands and a note is written. The rule never flips on that evidence.
3. CONFIRMATION by profile. The alpha-box / HMG-box gene-identity profiles
   (`db/Ascomycota/classifiers/MAT_core`) score each side's own protein on its own
   profile. If O's gene beats W's gene by `PROFILE_MIN_MARGIN` bits the call is
   O. If the profiles favour W, are inconclusive, or are unavailable, the call is
   `undetermined`. Neither the coverage rule nor the profiles decide alone.

Every triggered locus carries `idiomorph_vote_note` in the report, with both
reasons. Only non-optional core genes are considered, so the rule is independent
of `--idiomorph-policy` (which governs `optional` genes).

Everything here is pure; the pipeline supplies the per-gene facts.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

POLICY_OFF = "off"
POLICY_COVERAGE_PROFILE = "coverage_profile"
PARTIAL_CORE_POLICIES = (POLICY_OFF, POLICY_COVERAGE_PROFILE)
DEFAULT_PARTIAL_CORE_POLICY = POLICY_OFF

#: Percent of the reference protein a modelled gene covers (see
#: `pipeline.model_reference_coverage`). Both PROVISIONAL; the choice is argued
#: in docs/notes/2026-10-17_partial-coverage-vote.md from the v0.6.0 distribution.
PARTIAL_COVERAGE = 80.0
FULL_COVERAGE = 90.0
#: A partial hit this close to a contig end (either side) is an assembly artefact.
EDGE_BP = 500
#: Bits by which O's own profile must beat W's. PROVISIONAL: borrowed from the
#: Mucoromycota classifier's `min_margin: 25` (db/Mucoromycota/order.yml); the
#: Ascomycota alpha/HMG profiles have not had their own margin calibrated.
PROFILE_MIN_MARGIN = 25.0

#: Where a phylum keeps the gene-identity profiles, relative to its db directory.
CORE_PROFILE_SUBDIR = "classifiers/MAT_core"

DECISION_FLIPPED = "flipped"
DECISION_UNDETERMINED = "undetermined"
DECISION_KEPT_EDGE = "kept_edge_artefact"


def check_policy(policy: str) -> str:
    if policy not in PARTIAL_CORE_POLICIES:
        raise ValueError(f"unknown partial-core policy {policy!r}; expected one of {PARTIAL_CORE_POLICIES}")
    return policy


@dataclass(frozen=True)
class CoreStat:
    """What the rule needs to know of one modelled non-optional core gene."""

    gene: str
    idiomorph: str
    coverage: float | None  # percent of reference, None when unmodelled
    bitscore: float | None
    edge_distance: int | None  # bp to the nearer contig end, None when unknown


@dataclass(frozen=True)
class VoteDecision:
    decision: str  # DECISION_*
    idiomorph: str | None  # the call when flipped; None otherwise
    note: dict


def evaluate(
    vote: list[dict],
    stats: list[CoreStat],
    profile_scores: Callable[[list[str]], dict[str, float] | None] | None,
    *,
    partial: float = PARTIAL_COVERAGE,
    full: float = FULL_COVERAGE,
    edge_bp: int = EDGE_BP,
    margin: float = PROFILE_MIN_MARGIN,
) -> VoteDecision | None:
    """None when the rule does not trigger. Otherwise the decision and its note.

    `vote` is `idiomorph_candidates` output (best first). `profile_scores(genes)`
    returns each gene's best bit score on its own idiomorph's profile, or None
    when profiles are unavailable (missing files, no pyhmmer, no protein).
    """
    if len(vote) < 2 or not vote[0]["score"] > vote[1]["score"]:
        return None  # nothing to override: no contest, or an exact tie (already `undetermined`)
    winner, other = vote[0]["idiomorph"], vote[1]["idiomorph"]
    w_stats = [s for s in stats if s.idiomorph == winner and s.coverage is not None]
    o_stats = [s for s in stats if s.idiomorph == other and s.coverage is not None]
    if not w_stats or not o_stats:
        return None
    if any(s.coverage >= full for s in w_stats):
        return None  # the winner has a full-length core gene of its own
    w = max(w_stats, key=lambda s: (s.bitscore or 0.0))
    if w.coverage >= partial:
        return None
    o_full = [s for s in o_stats if s.coverage >= full]
    if not o_full:
        return None
    o = max(o_full, key=lambda s: (s.bitscore or 0.0))
    note: dict = {
        "rule": POLICY_COVERAGE_PROFILE,
        "winner_idiomorph": winner, "winner_gene": w.gene,
        "winner_coverage": w.coverage, "winner_bitscore": w.bitscore,
        "other_idiomorph": other, "other_gene": o.gene,
        "other_coverage": o.coverage, "other_bitscore": o.bitscore,
        "thresholds": {"partial_coverage": partial, "full_coverage": full,
                       "edge_bp": edge_bp, "profile_min_margin_bits": margin,
                       "profile_margin_provisional": True},
    }
    if w.edge_distance is not None and w.edge_distance < edge_bp:
        note.update(decision=DECISION_KEPT_EDGE, winner_edge_distance=w.edge_distance,
                    reason=(f"{w.gene} covers {w.coverage:.0f}% of its reference but lies {w.edge_distance} bp "
                            f"from a contig end: treated as truncation by the assembly, vote kept ({winner})"))
        return VoteDecision(DECISION_KEPT_EDGE, None, note)
    scores = profile_scores([w.gene, o.gene]) if profile_scores is not None else None
    if not scores or w.gene not in scores or o.gene not in scores:
        note.update(decision=DECISION_UNDETERMINED, profile=None,
                    reason=(f"coverage rule: {w.gene} is partial ({w.coverage:.0f}%) and {o.gene} is full length "
                            f"({o.coverage:.0f}%) -> {other}; profiles unavailable, so not confirmed"))
        return VoteDecision(DECISION_UNDETERMINED, None, note)
    diff = scores[o.gene] - scores[w.gene]
    profile_call = other if diff >= margin else winner if diff <= -margin else "inconclusive"
    note["profile"] = {"scores": {g: round(v, 1) for g, v in scores.items()},
                       "margin_bits": round(diff, 1), "call": profile_call}
    if profile_call == other:
        note.update(decision=DECISION_FLIPPED,
                    reason=(f"coverage rule ({w.gene} {w.coverage:.0f}% vs {o.gene} {o.coverage:.0f}% full length) "
                            f"and profiles ({o.gene} {scores[o.gene]:.0f} vs {w.gene} {scores[w.gene]:.0f} bits) "
                            f"agree: {other}, not the bitscore vote's {winner}"))
        return VoteDecision(DECISION_FLIPPED, other, note)
    note.update(decision=DECISION_UNDETERMINED,
                reason=(f"coverage rule favours {other} ({w.gene} {w.coverage:.0f}% partial, {o.gene} "
                        f"{o.coverage:.0f}% full) but profiles {'favour ' + winner if profile_call == winner else 'are inconclusive'} "
                        f"({o.gene} {scores[o.gene]:.0f} vs {w.gene} {scores[w.gene]:.0f} bits, margin "
                        f"{diff:+.0f}, needs {margin:g}); disagreement, so undetermined"))
    return VoteDecision(DECISION_UNDETERMINED, None, note)


def load_core_profiles(db_root: Path, family):
    """An `IdiomorphClassifier` over the family's gene-identity profiles, or None.

    The Ascomycota MAT_core directory is deliberately not named by `order.yml`
    (its profiles identify the GENE, alpha box or HMG box, and are not allele
    profiles), so the `idiomorph_classifier` machinery does not run for it. This
    loads the same files for the confirmation step only. None when the files, or
    pyhmmer, are absent.
    """
    from MATPredict.detect import classifier as clf_mod

    directory = Path(db_root) / family.key.phylum / CORE_PROFILE_SUBDIR
    if not (directory / "manifest.yaml").exists():
        return None
    try:
        import pyhmmer
        import hashlib
        models = {}
        for gene, idiomorph in clf_mod.classifier_genes(family).items():
            path = directory / f"{gene}.hmm"
            if path.exists():
                with pyhmmer.plan7.HMMFile(str(path)) as handle:
                    models[idiomorph] = (gene, handle.read())
        if len(models) < 2:
            return None
        return clf_mod.IdiomorphClassifier(
            directory=directory, models=models, min_margin=PROFILE_MIN_MARGIN,
            manifest_sha256=hashlib.sha256((directory / "manifest.yaml").read_bytes()).hexdigest())
    except ImportError:
        return None
