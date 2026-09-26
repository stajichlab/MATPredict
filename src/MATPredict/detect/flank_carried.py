"""The flank-carried rule: a call whose core genes were never modelled.

A flank-carried call cleared the modelled-gene bar on its FLANKING gene models
alone; every core_MAT gene in it is a bare alignment. Curator's ruling
2026-09-26 (docs/notes/2026-09-26_cauris-flank-carried-and-calbicans-
zygosity.md), measured in the Serinales-wide scan: 131 of 2,647 calls were
flank-carried, and 52 of them were Debaryomyces artefacts -- in CBS767 the
PAP1-OBP1-PIK1 block sits ~700 kb from the real a1/a2/alpha1 locus, and a 29%
MTLA2 fragment 8 kb from the block made a spurious second call.

Revised the same day, after the Ascomycota audit
(results/2026-09-26_flank_rule_ascomycota/NOTE.md), curator's ruling. The first
version required EVERY core hit inside the flank span +-3 kb. On the cap-6
panel it changed 89 calls, and it had two faults: one stray core hit 40-60 kb
away withheld a good call (Didymobotryum rigidum, whose MAT1-1-2, MAT1-1-3 and
MAT1-2-1 all sit between SLA2 and APN2), and 3 kb only fits families whose
flanks sit INSIDE the idiomorph (Serinales MTL); SLA2/APN2/COX13 sit outside
it. The rule now judges the call on its STRONGEST core hit (lowest e-value):

1. That hit has E <= `FLANK_CARRIED_MAX_EVALUE` and lies within the family's
   `flank_carried_window_bp` of the flank span, on the flank contig: the call
   is kept, capped at `low`, classed `partial_locus` and flagged
   `idiomorph_unmodelled`.
2. Otherwise (a weak or noise-level strongest hit, too far, another contig, or
   no core hit at all): the call is withheld, exactly like a bar failure.

Simulated on the audit's 108 changed calls: Ascomycota keeps 9 of 89 at low,
the Mucoromycota group 8 of 19, and all 3 real loci the first version
withheld (two Mucor irregularis, Trigonopsis variabilis) are kept.

Generic over each family's own flanking genes (any `flanking_*` role), not
hardcoded to the Serinales PAP1/OBP1/PIK1 roster. The span is taken over EVERY
flank hit in the call, modelled or not, matching `results/2026-09-26_
serinales_all_882aa01/flank_carried_audit.py`.
"""
from __future__ import annotations

import math
from dataclasses import replace
from typing import Mapping

from MATPredict.detect.family_registry import DEFAULT_FLANK_CARRIED_WINDOW_BP
from MATPredict.detect.idiomorph import LOCUS_CLASS_PARTIAL
from MATPredict.detect.polish import STATUS_AGREE, STATUS_DISAGREE, STATUS_SINGLE

#: The strongest core hit must reach this e-value for a flank-carried call to
#: be kept. NOT tuned: it sits in the wide gap the audit measured between
#: strong hits (E <= 1e-5) and noise (E > 1e-2) on the 108 changed calls.
FLANK_CARRIED_MAX_EVALUE = 1e-5

#: `DetectionResult.withheld_reason` for a call this rule withholds.
WITHHELD_FLANK_CARRIED = "flank_carried_core_outside_flank_span"

_MODELLED_STATUSES = (STATUS_AGREE, STATUS_DISAGREE, STATUS_SINGLE)


def _modelled(evidence) -> bool:
    """A real gene model: polished by a tool, or annotated (the fast path).
    The same two sources `_modelled_gene_names` counts toward the bar."""
    return evidence.status in _MODELLED_STATUSES or evidence.method == "diamond_proteome"


def _is_flank(evidence) -> bool:
    return evidence.role.startswith("flanking")


def _evalue(evidence) -> float:
    """A missing e-value ranks as the weakest possible hit, never the strongest."""
    value = getattr(evidence, "evalue", None)
    return math.inf if value is None else value


def apply_flank_carried_rule(
    results: list,
    window_bp_by_family: Mapping | None = None,
    max_evalue: float = FLANK_CARRIED_MAX_EVALUE,
):
    """Split `results` into (kept, withheld) under the flank-carried rule.

    `window_bp_by_family` maps a `FamilyKey` to its `flank_carried_window_bp`;
    a family missing from it gets `DEFAULT_FLANK_CARRIED_WINDOW_BP`.

    A result with a modelled core gene, or with no modelled flank, is not
    flank-carried and passes through unchanged: the modelled-gene bar and
    tiering already rule on it.
    """
    windows = window_bp_by_family or {}
    kept, withheld = [], []
    for result in results:
        core = [e for e in result.gene_evidence if e.role == "core_MAT"]
        flanks = [e for e in result.gene_evidence if _is_flank(e)]
        if any(_modelled(e) for e in core) or not any(_modelled(e) for e in flanks):
            kept.append(result)
            continue
        contig = flanks[0].contig
        low = min(e.start for e in flanks if e.contig == contig)
        high = max(e.end for e in flanks if e.contig == contig)

        def distance(e) -> float:
            if e.contig != contig:
                return math.inf
            return max(low - e.end, e.start - high, 0)

        window = windows.get(result.family_key, DEFAULT_FLANK_CARRIED_WINDOW_BP)
        # Strongest = lowest e-value; the nearer hit breaks a tie.
        strongest = min(core, key=lambda e: (_evalue(e), distance(e))) if core else None
        if (strongest is not None and _evalue(strongest) <= max_evalue
                and distance(strongest) <= window):
            kept.append(replace(
                result, confidence="low", locus_class=LOCUS_CLASS_PARTIAL,
                idiomorph_unmodelled=True,
            ))
        else:
            withheld.append(replace(result, withheld_reason=WITHHELD_FLANK_CARRIED))
    return kept, withheld
