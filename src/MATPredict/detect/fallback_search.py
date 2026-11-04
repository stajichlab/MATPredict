"""Search the phylum-fallback families when a routed family's hits were all withheld (#228).

A genome routed by taxid or lineage searches only the families whose scope contains it (Ogataea: `MATog`). If that search reports no
locus and the modelled-gene bar withheld at least one, v0.7.0 would have called the genome from the phylum-fallback families
(`MATsc`, `MTL`); v0.8.0 leaves it uncalled. With the policy `on` those families are searched as well, only in that case, and their
calls go through the usual rules for a family outside its scope (`family_status`). The first search is not repeated. `off` is the
default, so the v0.7 and v0.8 presets reproduce their releases.
"""
from __future__ import annotations

from dataclasses import replace

from MATPredict.detect.pipeline import MODELLED_GENE_BAR, DetectionOutcome

POLICY_OFF = "off"
POLICY_ON = "on"
FALLBACK_AFTER_WITHHELD_POLICIES = (POLICY_OFF, POLICY_ON)
DEFAULT_FALLBACK_AFTER_WITHHELD_POLICY = POLICY_OFF

#: Routing modes of a search that names its families by scope; the others (`phylum_fallback`, `explicit_phylum`, `exhaustive`) already
#: search the phylum, so there is nothing to add.
SCOPED_ROUTES = ("direct", "lineage")


def needed(outcome: DetectionOutcome) -> bool:
    """True when the routed search reported nothing and the modelled-gene bar withheld something."""
    return (outcome.routing_mode in SCOPED_ROUTES and not outcome.results
            and any(r.withheld_reason == MODELLED_GENE_BAR for r in outcome.suppressed_loci))


def merge(first: DetectionOutcome, second: DetectionOutcome, second_routing_mode: str) -> DetectionOutcome:
    """`first` plus what the fallback families found; `routing_mode` stays the routed one and `fallback_search` says what was added."""
    return replace(
        first,
        results=list(first.results) + list(second.results),
        not_detected=list(first.not_detected) + list(second.not_detected),
        families_attempted=list(first.families_attempted) + [k for k in second.families_attempted if k not in first.families_attempted],
        suppressed_unpolished=first.suppressed_unpolished + second.suppressed_unpolished,
        suppressed_loci=list(first.suppressed_loci) + list(second.suppressed_loci),
        suppressed_flank_carried=first.suppressed_flank_carried + second.suppressed_flank_carried,
        suppressed_mat_gene_gate=first.suppressed_mat_gene_gate + second.suppressed_mat_gene_gate,
        suppressed_paralog_class=first.suppressed_paralog_class + second.suppressed_paralog_class,
        assembly_gaps_at_locus=list(first.assembly_gaps_at_locus) + list(second.assembly_gaps_at_locus),
        two_idiomorphs=list(first.two_idiomorphs) + list(second.two_idiomorphs),
        fallback_search={
            "reason": "routed family reported nothing: every hit was withheld by the modelled-gene bar",
            "routing_mode": second_routing_mode,
            "families": [f"{k.phylum}:{k.locus_name}" for k in second.families_attempted],
            "calls": len(second.results),
        },
    )
