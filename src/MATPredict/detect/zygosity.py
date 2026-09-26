"""Genome-level zygosity for taxa whose assemblies collapse MTL heterozygosity.

Curator's ruling 2026-09-26 (results/2026-09-26_calbicans_mtl_reads/NOTE.md).
Read depth over the C. albicans MTL idiomorphs split 16 isolates cleanly into
heterozygotes and homozygotes. 8 of the 11 read-based a/alpha isolates had an
assembly that holds only one idiomorph; no assembly showed both when the reads
showed one. So in such a taxon an assembly-based single-idiomorph call cannot
tell a homozygote from a collapsed heterozygote: zygosity is unknown.

The taxa are curated data in `db/assembly_zygosity.yml`, each with its reason
and the evidence behind it, so adding a species is a data change. The rule
never alters an idiomorph call; it only adds a genome-level statement.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import yaml

#: The curated list, relative to the database root.
ZYGOSITY_FILE = "assembly_zygosity.yml"

_NO_IDIOMORPH = {None, "undetermined"}


@dataclass(frozen=True)
class ZygosityRule:
    taxid: int
    name: str
    reason: str
    evidence: str


def load_zygosity_rules(db_root: Path) -> list[ZygosityRule]:
    """The curated taxa, or [] when the file does not exist."""
    path = Path(db_root) / ZYGOSITY_FILE
    if not path.exists():
        return []
    doc = yaml.safe_load(path.read_text()) or {}
    return [
        ZygosityRule(int(e["taxid"]), e["name"], e["reason"], e["evidence"])
        for e in doc.get("assemblies_collapse_heterozygosity") or []
    ]


def genome_zygosity(
    results: list, taxid: int | None, rules: list[ZygosityRule],
    lineage_resolver: Callable[[int], list[int]],
) -> dict | None:
    """The genome's zygosity statement, or None when no rule applies.

    Applies when the genome's reported calls name exactly one idiomorph and
    its taxid, or an ancestor, is a listed taxon. The lineage is looked up
    only when that is the last open question, and a failed lookup gives
    `unchecked` with the error rather than a silent pass.
    """
    idiomorphs = {r.idiomorph for r in results} - _NO_IDIOMORPH
    if taxid is None or not rules or len(idiomorphs) != 1:
        return None
    by_taxid = {r.taxid: r for r in rules}
    rule = by_taxid.get(taxid)
    if rule is None:
        try:
            lineage = lineage_resolver(taxid)
        except Exception as exc:  # noqa: BLE001 - recorded in the report, never fatal
            return {
                "status": "unchecked",
                "reason": f"lineage lookup failed ({type(exc).__name__}: {str(exc)[:160]}); "
                          "cannot tell whether this taxon's assemblies collapse "
                          "MTL heterozygosity",
                "taxid": None,
                "evidence": None,
            }
        rule = next((by_taxid[t] for t in lineage if t in by_taxid), None)
    if rule is None:
        return None
    return {
        "status": "unknown",
        "reason": f"{rule.name}: {rule.reason}; a single-idiomorph assembly call "
                  "cannot distinguish a homozygote from a collapsed heterozygote",
        "taxid": rule.taxid,
        "evidence": rule.evidence,
    }
