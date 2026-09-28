"""Withhold a second call the classifier could not assign.

Curator's ruling 2026-09-27. Adding the tier-2 Umbelopsis records let weak
flank hits (28-38% identity) be modelled next to HMG paralogs, and 13
Mucoromycota genomes gained an 85-155 kb second call whose idiomorph the HMM
classifier scored but could not assign (margin 19.5-24.8, below the family's
25-bit `min_margin`), while the genome already had a determined call of the
same family (results/2026-09-27_umbelopsis_curation/;
results/2026-09-27_mucoro_curation_guard/replay_guard.txt).

Rule: a reported call is withheld when ALL hold --

1. the idiomorph classifier ran on it (`idiomorph_classifier` is set) and
   returned `undetermined`;
2. another reported call of the SAME family in the genome has a determined
   idiomorph;
3. it is not a `homothallic_candidate` (such a call rests on both alleles).

Narrow on purpose. A replay of a broader version (every enum-vocabulary
family) would have withheld 8 Saccharomyces silent-cassette calls
(MATALPHA2/MATA2 only), where several loci per genome are normal; those have
no classifier verdict, so condition 1 leaves them alone. A lone undetermined
call (nothing to defer to) is kept.
"""
from __future__ import annotations

from dataclasses import replace

#: `DetectionResult.withheld_reason` for a call this rule withholds.
WITHHELD_SECONDARY_UNDETERMINED = "secondary_call_classifier_undetermined"

_HOMOTHALLIC = "homothallic_candidate"


def _determined(result) -> bool:
    return result.idiomorph not in (None, "undetermined")


def _classifier_undetermined(result) -> bool:
    verdict = result.idiomorph_classifier
    return (
        verdict is not None
        and result.idiomorph == "undetermined"
        and verdict.get("idiomorph") == "undetermined"
    )


def apply_secondary_undetermined_rule(results):
    """Split `results` into (kept, withheld) under the rule above."""
    families_determined = {r.family_key for r in results if _determined(r)}
    kept, withheld = [], []
    for r in results:
        if (
            _classifier_undetermined(r)
            and r.locus_class != _HOMOTHALLIC
            and r.family_key in families_determined
        ):
            withheld.append(replace(r, withheld_reason=WITHHELD_SECONDARY_UNDETERMINED))
        else:
            kept.append(r)
    return kept, withheld
