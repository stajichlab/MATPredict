"""Label a call made by searching a family outside its curated phylum.

Curator's ruling 2026-09-26: the Mortierellomycota and Kickxellomycota calls,
made with `--phylum Mucoromycota`, are "unverified -- not clear what they are
yet". The flank-ortholog synteny check (results/2026-09-26_flank_synteny_ED/
NOTE.md) found the Mucorales tptA-HMG-rnhA arrangement in 0/100
Mortierellomycota and 0/190 Kickxellomycota genomes, against 117/293 in the
Mucoromycota control: none of their 13 calls is supported.

Only an override route can search a family outside the genome's phylum:
`explicit_phylum` (`--phylum`) or `exhaustive` (`--exhaustive`). Taxid
routing searches the genome's own phylum by construction and is never
labelled. On an override route a call is unverified when the genome's phylum
is known and differs from the family's, or -- on `exhaustive` only -- when the
genome's phylum is unknown, since nothing then ties the genome to any curated
phylum. An `explicit_phylum` run with an unknown genome phylum is not
labelled: the operator stated the phylum and nothing contradicts it.

The label never changes a call's confidence, class or idiomorph.
"""
from __future__ import annotations

import logging
from dataclasses import replace


logger = logging.getLogger("MATPredict")

#: Routes on which a family can be searched outside the genome's phylum.
OVERRIDE_ROUTES = ("explicit_phylum", "exhaustive")

#: The measurement the label rests on.
UNVERIFIED_EVIDENCE = "results/2026-09-26_flank_synteny_ED/NOTE.md"


def label_verification(results: list, routing_mode: str | None, genome_phylum: str | None):
    """`results` with `verification` set on every call the ruling labels."""
    if routing_mode not in OVERRIDE_ROUTES:
        return list(results)
    out = []
    for r in results:
        family_phylum = r.family_key.phylum
        if genome_phylum is not None and genome_phylum == family_phylum:
            out.append(r)
            continue
        if genome_phylum is None and routing_mode == "explicit_phylum":
            out.append(r)
            continue
        where = (f"a {genome_phylum} genome" if genome_phylum
                 else "a genome of unknown phylum")
        out.append(replace(r, verification={
            "status": "unverified",
            "reason": (
                f"{family_phylum} family searched on {where} by the "
                f"{routing_mode} override; the curated locus architecture is "
                "not shown to hold outside its phylum"
            ),
            "family_phylum": family_phylum,
            "genome_phylum": genome_phylum,
            "evidence": UNVERIFIED_EVIDENCE,
        }))
    return out


# ---------------------------------------------------------------------------
# CAAX-dependent calls (curator's ruling 2026-09-28, review finding F4; it
# replaces the 2026-09-27 four-family list). A call that reaches the admission
# bar only by counting a strict-CAAX scan precursor
# (`DetectionResult.caax_dependent`) is unverified wherever it occurs: the
# negative control could not bound the finder's false-positive rate (6/9
# curated mating receptors flagged, 0/25 non-mating STE3 copies with a 95%
# upper limit of 13.7%, random windows 2.3%). Never changes confidence.
# REVIEW LATER, once a labelled set of >= 100 non-mating STE3 loci exists.

#: The measurement behind the label, relative to the repository root.
CAAX_UNVERIFIED_EVIDENCE = "results/2026-09-28_validation_f3_f4/NOTE.md"


def label_caax_unverified(results: list) -> list:
    """`results` with every CAAX-dependent call labelled `unverified`. A call
    that already carries a `verification` label keeps it."""
    label = {
        "status": "unverified",
        "reason": (
            "admitted only through a strict-CAAX scan precursor; the scan's "
            "false-positive rate is not yet bounded (0/25 non-mating STE3 copies "
            "flagged, 95% upper limit 13.7%)"
        ),
        "evidence": CAAX_UNVERIFIED_EVIDENCE,
    }
    return [replace(r, verification=dict(label)) if (r.caax_dependent and r.verification is None)
            else r for r in results]
