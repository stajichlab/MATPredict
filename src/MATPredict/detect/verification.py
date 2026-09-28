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
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable

import yaml

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
# CAAX-dependent PR calls in low-enrichment families (curator's ruling
# 2026-09-27). A Basidiomycota:PR call that reaches the admission bar only by
# counting a strict-CAAX scan precursor (`DetectionResult.caax_dependent`) is
# unverified when the genome lies in a curated taxon where the CAAX-gained
# calls showed low enrichment for the curated mating-receptor clades
# (results/2026-09-27_caax_precursor/per_family.tsv, gained.tsv). The taxa are
# curation data in `db/caax_unverified_taxa.yml`. Never changes confidence.


#: The curated taxa, relative to the database root.
CAAX_UNVERIFIED_FILE = "caax_unverified_taxa.yml"


@dataclass(frozen=True)
class CaaxUnverifiedRule:
    taxid: int
    name: str
    reason: str
    evidence: str


def load_caax_unverified_rules(db_root: Path) -> list[CaaxUnverifiedRule]:
    """The curated taxa, or [] when the file does not exist."""
    path = Path(db_root) / CAAX_UNVERIFIED_FILE
    if not path.exists():
        return []
    doc = yaml.safe_load(path.read_text()) or {}
    return [
        CaaxUnverifiedRule(int(e["taxid"]), e["name"], e["reason"], e["evidence"])
        for e in doc.get("caax_unverified_taxa") or []
    ]


def label_caax_unverified(
    results: list, taxid: int | None, rules: list[CaaxUnverifiedRule],
    lineage_resolver: Callable[[int], list[int]],
) -> list:
    """`results` with CAAX-dependent calls in a listed taxon labelled.

    The genome matches when its taxid, or any ancestor, is a listed taxid.
    The lineage is looked up only when a CAAX-dependent call exists and the
    genome's own taxid is not listed; a failed lookup leaves every call
    unlabelled (logged), never an error. A call that already carries a
    `verification` label keeps it.
    """
    targets = [r for r in results if r.caax_dependent and r.verification is None]
    if taxid is None or not rules or not targets:
        return list(results)
    by_taxid = {r.taxid: r for r in rules}
    rule = by_taxid.get(taxid)
    if rule is None:
        try:
            lineage = lineage_resolver(taxid) or []
        except Exception as exc:  # noqa: BLE001 - logged, no label rather than a lost call
            logger.warning("lineage lookup for the CAAX unverified label failed: %s", exc)
            return list(results)
        rule = next((by_taxid[t] for t in lineage if t in by_taxid), None)
    if rule is None:
        return list(results)
    label = {
        "status": "unverified",
        "reason": (
            f"admitted only through a strict-CAAX scan precursor in {rule.name}, "
            f"where {rule.reason}"
        ),
        "taxid": rule.taxid,
        "taxon": rule.name,
        "evidence": rule.evidence,
    }
    return [replace(r, verification=label) if (r.caax_dependent and r.verification is None) else r
            for r in results]
