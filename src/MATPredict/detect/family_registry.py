"""Load order.yml families and route them to a taxid via taxonomic_scope."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import yaml

from MATPredict.db.taxonomy import TaxonomyResult, resolve_lineage


@dataclass(frozen=True)
class FamilyKey:
    phylum: str
    locus_name: str


@dataclass(frozen=True)
class Family:
    key: FamilyKey
    vocabulary_type: str
    idiomorph_values: list[str] | None
    idiomorph_pattern: str | None
    genes: list[dict]
    taxonomic_scope: list[int]


def load_all_families(db_root: Path) -> list[Family]:
    """Read every db/<Phylum>/order.yml and flatten it into Family records."""
    families: list[Family] = []
    for order_file in sorted(db_root.glob("*/order.yml")):
        doc = yaml.safe_load(order_file.read_text())
        for locus in doc["loci"]:
            families.append(
                Family(
                    key=FamilyKey(doc["phylum"], locus["locus_name"]),
                    vocabulary_type=locus["vocabulary_type"],
                    idiomorph_values=locus.get("idiomorph_values"),
                    idiomorph_pattern=locus.get("idiomorph_pattern"),
                    genes=locus["genes"],
                    taxonomic_scope=locus["taxonomic_scope"],
                )
            )
    return families


def _lineage_contains_taxid(lineage: str, taxid: int) -> bool:
    # taxonkit's k__/p__/... string doesn't carry per-rank taxids, so an exact
    # numeric match against the scope list is only possible for the queried
    # taxid itself; broader ancestor matching needs the full taxid lineage,
    # not the name-string lineage taxonomy.resolve_lineage returns today.
    # For v1, scope membership is: the queried taxid is itself in scope, OR
    # a scope taxid's name appears in the lineage string (best-effort).
    return False  # overwritten by route() below using the numeric taxid directly


def route(
    taxid: int | None,
    families: list[Family],
    lineage_resolver: Callable[[int], TaxonomyResult] = resolve_lineage,
) -> list[Family]:
    """Return families whose taxonomic_scope contains taxid.

    v1 scope matching is exact-taxid-membership only (a family's
    taxonomic_scope must directly list the queried taxid, since
    MATPredict.db.taxonomy.resolve_lineage returns a rank-name string, not a
    numeric ancestor chain, so ancestor-subtree matching isn't available yet).
    taxid=None, or a taxid matching no family's scope, returns every family
    unchanged (the exhaustive fallback path)."""
    if taxid is None:
        return list(families)
    lineage_resolver(taxid)  # resolved for future ancestor-aware matching; unused in v1 matching itself
    matched = [f for f in families if taxid in f.taxonomic_scope]
    return matched if matched else list(families)
