"""The flank-carried rule: a call whose core genes were never modelled.

A flank-carried call cleared the modelled-gene bar on its FLANKING gene models
alone; every core_MAT gene in it is a bare alignment. Curator's ruling
2026-09-26 (docs/notes/2026-09-26_cauris-flank-carried-and-calbicans-
zygosity.md), measured in the Serinales-wide scan: 131 of 2,647 calls were
flank-carried, and 52 of them were Debaryomyces artefacts -- in CBS767 the
PAP1-OBP1-PIK1 block sits ~700 kb from the real a1/a2/alpha1 locus, and a 29%
MTLA2 fragment 8 kb from the block made a spurious second call.

1. Every core hit lies inside the flank span, padded by
   `FLANK_SPAN_PADDING_BP` on each side: the call is kept, capped at `low`,
   classed `partial_locus` and flagged `idiomorph_unmodelled` -- a real locus
   position whose idiomorph rests on unmodelled hits.
2. Any core hit lies outside it (or on another contig, or there is no core
   hit at all): the call is withheld, exactly like a bar failure.

Generic over each family's own flanking genes (any `flanking_*` role), not
hardcoded to the Serinales PAP1/OBP1/PIK1 roster, so it applies wherever a
family curates flanks. The span is taken over EVERY flank hit in the call,
modelled or not, matching `results/2026-09-26_serinales_all_882aa01/
flank_carried_audit.py`, which produced the measurement the ruling rests on.
"""
from __future__ import annotations

from dataclasses import replace

from MATPredict.detect.idiomorph import LOCUS_CLASS_PARTIAL
from MATPredict.detect.polish import STATUS_AGREE, STATUS_DISAGREE, STATUS_SINGLE

#: Padding on each side of the flank span, from the curator's ruling (+-3 kb).
FLANK_SPAN_PADDING_BP = 3_000

#: `DetectionResult.withheld_reason` for a call this rule withholds.
WITHHELD_FLANK_CARRIED = "flank_carried_core_outside_flank_span"

_MODELLED_STATUSES = (STATUS_AGREE, STATUS_DISAGREE, STATUS_SINGLE)


def _modelled(evidence) -> bool:
    """A real gene model: polished by a tool, or annotated (the fast path).
    The same two sources `_modelled_gene_names` counts toward the bar."""
    return evidence.status in _MODELLED_STATUSES or evidence.method == "diamond_proteome"


def _is_flank(evidence) -> bool:
    return evidence.role.startswith("flanking")


def apply_flank_carried_rule(results: list, padding_bp: int = FLANK_SPAN_PADDING_BP):
    """Split `results` into (kept, withheld) under the flank-carried rule.

    A result with a modelled core gene, or with no modelled flank, is not
    flank-carried and passes through unchanged: the modelled-gene bar and
    tiering already rule on it.
    """
    kept, withheld = [], []
    for result in results:
        core = [e for e in result.gene_evidence if e.role == "core_MAT"]
        flanks = [e for e in result.gene_evidence if _is_flank(e)]
        if any(_modelled(e) for e in core) or not any(_modelled(e) for e in flanks):
            kept.append(result)
            continue
        contig = flanks[0].contig
        low = min(e.start for e in flanks if e.contig == contig) - padding_bp
        high = max(e.end for e in flanks if e.contig == contig) + padding_bp
        inside = bool(core) and all(
            e.contig == contig and low <= e.start and e.end <= high for e in core
        )
        if inside:
            kept.append(replace(
                result, confidence="low", locus_class=LOCUS_CLASS_PARTIAL,
                idiomorph_unmodelled=True,
            ))
        else:
            withheld.append(replace(result, withheld_reason=WITHHELD_FLANK_CARRIED))
    return kept, withheld
