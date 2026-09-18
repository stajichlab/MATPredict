"""Load order.yml families and route them to a taxid via taxonomic_scope."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import yaml

from MATPredict.db.taxonomy import default_lineage_taxids


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


def load_record_families(db_root: Path) -> dict[str, FamilyKey]:
    """Map every accepted curated record_id to the one family it belongs to.

    A curated record lives at `db/<Phylum>/<Order-or-Family>/<record_id>/metadata.yaml`
    and declares exactly one `mating_type.locus_name`, so `(phylum, locus_name)` --
    i.e. its `FamilyKey` -- is unambiguous per record. This index is what lets
    `search.py` attribute a hit to the family whose curated protein it actually
    matched, instead of guessing from the bare gene name (gene names such as
    `pheromone`, `pheromone_receptor`, `Z`, `Y`, `matPc`, `matMc` and `sla2` are
    reused by several distinct families in the real database, so a bare-name
    lookup silently collapses those families into whichever one is written last).

    `db/candidates/...` matches the same glob shape but holds proposed, not
    accepted, records; it is excluded by name exactly as `benchmark._load_records`
    does.
    """
    index: dict[str, FamilyKey] = {}
    for meta_path in sorted(db_root.glob("*/*/*/metadata.yaml")):
        if meta_path.relative_to(db_root).parts[0] == "candidates":
            continue
        doc = yaml.safe_load(meta_path.read_text())
        record_id = doc.get("record_id")
        locus_name = (doc.get("mating_type") or {}).get("locus_name")
        if not record_id or not locus_name:
            continue
        index[record_id] = FamilyKey(meta_path.parents[2].name, locus_name)
    return index


def route(
    taxid: int | None,
    families: list[Family],
    lineage_taxids_resolver: Callable[[int], list[int]] = default_lineage_taxids,
) -> list[Family]:
    """Return families whose taxonomic_scope contains taxid, directly or via lineage.

    A family matches if EITHER the queried taxid is directly listed in its
    taxonomic_scope (the original exact-membership check) OR the taxid's NCBI
    Taxonomy ancestor lineage contains any taxid in its taxonomic_scope. Most
    families declare a broad scope (e.g. a subphylum/subclass taxid) expecting
    it to cover every descendant species -- lineage matching is what actually
    makes that work; before this, only records whose scope also happened to
    list their exact species/strain taxid routed correctly (an audit found 50
    of 61 curated records fell through to the exhaustive fallback below).

    The direct-membership check runs first and short-circuits before any
    lineage lookup (no network/subprocess call) whenever it already finds a
    match, so this stays free for the common case. `lineage_taxids_resolver`
    defaults to `MATPredict.db.taxonomy.default_lineage_taxids`, which fetches
    the ancestor chain via a cached NCBI Taxonomy efetch call; a resolver
    failure (network error, unknown taxid, etc.) degrades gracefully to the
    exhaustive fallback -- same as taxid=None or no scope matching at all.
    """
    if taxid is None:
        return list(families)
    direct = [f for f in families if taxid in f.taxonomic_scope]
    if direct:
        return direct
    try:
        ancestors = set(lineage_taxids_resolver(taxid))
    except Exception:
        ancestors = set()
    lineage_matched = [f for f in families if ancestors & set(f.taxonomic_scope)]
    return lineage_matched if lineage_matched else list(families)
