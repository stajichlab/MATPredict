"""Load order.yml families and route them to a taxid via taxonomic_scope."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import yaml

from MATPredict import logger
from MATPredict.db.taxonomy import default_lineage_phylum_name, default_lineage_taxids


DEFAULT_MAX_CLUSTER_GAP_BP = 25_000
"""Cluster gap used by any locus whose `order.yml` entry does not declare one.

25 kb was the single global constant this pipeline used for every phylum. The
curator ruled on 2026-09-20 that it is right for Ascomycota, so it stays the
default rather than becoming a value every locus must restate: absence in
`order.yml` means "the Ascomycota-calibrated 25 kb is fine here". Only a locus
that is known to need something else carries `max_cluster_gap_bp`.
"""


DEFAULT_MIN_IDIOMORPH_MARGIN = 5.0
"""Identity points two idiomorphs must be apart before the call is trusted.

Below this, the locus is still called -- refusing would cost real detections --
but its confidence tier is capped, and the margin is written into the report so
a close call is never mistaken for a clean one.

PROVISIONAL. Measured on the 23 ground-truth Mucoromycota genomes, where the
identity rule is 23/23 correct but the margins split hard by idiomorph: Plus
calls separate by 44.6-64.2 points, Minus calls by only 2.3-13.0. The narrowest
call seen anywhere (Blakeslea trispora, 0.59) sits outside the truth set, so
whether it is right is unknown. 5.0 caps the four calls under five points and
leaves the well-separated ones alone. Revise it from the resolution events in
the evidence diagnostics, not by intuition.
"""


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
    max_cluster_gap_bp: int = DEFAULT_MAX_CLUSTER_GAP_BP
    """How far apart two hits of this family may be and still be one locus.

    This is curation data, not a tuning knob, which is why it lives on the
    family (i.e. on the `order.yml` locus entry) instead of staying the single
    global constant it used to be. How spread out a MAT locus is, is a property
    of the clade's locus architecture, and it differs between clades: the
    curator ruled on 2026-09-20 that 25 kb is right for Ascomycota but too
    tight for Mucoromycota. Defaulted here so every existing construction --
    and every locus that does not declare one -- keeps the 25 kb behaviour
    exactly.
    """
    min_idiomorph_margin: float = DEFAULT_MIN_IDIOMORPH_MARGIN
    """Identity points two mutually exclusive idiomorph genes must be apart.

    Curation data on the same footing as `max_cluster_gap_bp`: how far apart
    two idiomorphs' identities fall depends on how densely that phylum's
    idiomorphs are represented in the curated database, and that differs
    sharply between phyla. A call below this margin is still made, but its
    tier is capped. See `DEFAULT_MIN_IDIOMORPH_MARGIN` for the measured basis.
    """


@dataclass(frozen=True)
class RoutingDecision:
    """Which families one detection run will search, and WHY they were chosen.

    `routing_mode` is returned rather than left for the caller to re-derive
    because the "why" changes how much the run's negatives are worth: an
    `exhaustive` run searched families from phyla the query cannot belong to,
    so its not-detected entries for those families are meaningless, whereas a
    `direct` run's are real evidence of absence. `report.write_detection_report`
    writes it into every detection report for exactly that reason.

    * `direct` -- the query taxid is itself listed in the matched families'
      `taxonomic_scope`.
    * `lineage` -- an ancestor of the query taxid is.
    * `phylum_fallback` -- nothing matched, so the run was narrowed to the
      query taxon's own phylum.
    * `exhaustive` -- nothing matched and the phylum could not be used, so
      every family in every phylum is searched. This is the expensive,
      low-value path; it is what Task 1 exists to make rare.
    * `explicit_phylum` -- the operator passed `--phylum`, overriding all
      taxid-based routing.

    `phylum` is set only for the two phylum-scoped modes, and is None otherwise.
    """

    families: list[Family]
    routing_mode: str
    phylum: str | None = None


def available_phyla(db_root: Path) -> list[str]:
    """The phylum names `db_root` actually contains, sorted, read at call time.

    Discovered from the `phylum:` field of each `db/<Phylum>/order.yml` -- the
    same files `load_all_families` reads and the same string that ends up in
    `FamilyKey.phylum` -- rather than from a hardcoded list, so adding a fourth
    phylum directory to the database makes it selectable with no code change,
    and so the `--phylum` CLI choices can never drift from what is on disk.
    Reading the field (not the directory name) guarantees the returned strings
    compare equal to `FamilyKey.phylum`; a directory with no `order.yml`
    (`db/candidates/`, `db/_schema/`) declares no phylum and is not offered.

    An `order.yml` that cannot be read or parsed costs its own phylum a
    `--phylum` choice and nothing else: it is logged by name and skipped, never
    raised. This function runs while the top-level argparse parser is being
    built, so raising here would abort `matpredict --help` and every subcommand
    that has nothing to do with the broken file -- a curator mid-edit would
    take down the whole CLI. Commands that genuinely need the file's contents
    (`load_all_families`) still fail loudly on it.
    """
    names = set()
    for order_file in db_root.glob("*/order.yml"):
        try:
            doc = yaml.safe_load(order_file.read_text())
        except (OSError, yaml.YAMLError) as err:
            logger.warning(
                "available_phyla: could not read %s (%s) -- that phylum will not be "
                "offered as a --phylum choice", order_file, err,
            )
            continue
        if doc and doc.get("phylum"):
            names.add(doc["phylum"])
    return sorted(names)


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
                    max_cluster_gap_bp=locus.get(
                        "max_cluster_gap_bp", DEFAULT_MAX_CLUSTER_GAP_BP
                    ),
                    min_idiomorph_margin=locus.get(
                        "min_idiomorph_margin", DEFAULT_MIN_IDIOMORPH_MARGIN
                    ),
                )
            )
    return families


def derive_max_cluster_gap(families: list[Family]) -> int:
    """The clustering gap for a run over `families`: the MAXIMUM of their gaps.

    The maximum, not the minimum or a per-family value, because the two errors
    are not symmetric. Taking too LARGE a gap under-splits -- two neighbouring
    loci can be merged into one cluster -- and that is recoverable downstream:
    the evidence floor and the polish stage still discriminate gene by gene
    within an over-large cluster, so the real locus is still there to be
    scored. Taking too SMALL a gap over-splits, cutting one real locus in two,
    and nothing downstream can put it back: each half is scored as an
    independent, incomplete candidate and the real locus is silently gone.

    A per-family gap is not possible here without changing `cluster_hits`,
    which groups hits by CONTIG ONLY and is deliberately family-agnostic
    (a real locus's hits are attributed to whichever curated family's protein
    they matched, so one locus's cluster routinely mixes families). One gap
    per run is therefore the unit of choice, and the maximum is the safe end.

    With no families (nothing routed) there is nothing to derive from, so the
    default stands.
    """
    return max((f.max_cluster_gap_bp for f in families), default=DEFAULT_MAX_CLUSTER_GAP_BP)


def expected_genes_for_idiomorph(
    family: Family, found_gene_names: set[str] | list[str]
) -> list[dict]:
    """The subset of `family.genes` actually expected given which idiomorph(s)
    the FOUND genes imply, so a real single-idiomorph genome is scored against
    only the genes that idiomorph should have -- not the family's full,
    multi-idiomorph gene roster.

    A gene with no `present_in_idiomorphs` (e.g. a flanking gene) applies to
    every idiomorph and is always included. When the found genes' own
    `present_in_idiomorphs` values name exactly ONE idiomorph, the returned
    list is narrowed to that idiomorph's genes plus every idiomorph-agnostic
    gene. When they name TWO OR MORE idiomorphs (a real, legitimate
    homothallic both-idiomorphs-present locus, already documented elsewhere
    in this project's curated records) or when none of the found genes carry
    a `present_in_idiomorphs` value at all (nothing to narrow by), the FULL
    roster is returned unchanged -- this is deliberately the conservative,
    already-existing behavior for those two cases, not a regression.
    """
    found = set(found_gene_names)
    found_idiomorphs: set[str] = set()
    for gene in family.genes:
        if gene["name"] in found:
            found_idiomorphs.update(gene.get("present_in_idiomorphs") or [])

    if len(found_idiomorphs) != 1:
        return family.genes  # ambiguous/both-present/nothing-to-narrow-by: full roster, unchanged behavior

    idiomorph = next(iter(found_idiomorphs))
    return [
        gene
        for gene in family.genes
        if not gene.get("present_in_idiomorphs") or idiomorph in gene["present_in_idiomorphs"]
    ]


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
    phylum_name_resolver: Callable[[int], str | None] = default_lineage_phylum_name,
    phylum: str | None = None,
) -> RoutingDecision:
    """Choose the families a detection run will search, and report which rule chose them.

    The rules are tried in this order, each one narrower and cheaper to search
    than the next:

    1. **`phylum` given** (`matpredict detect --phylum <P>`): an outright
       operator restriction to that phylum's families. No direct check, no
       lineage fetch, no fallback -- the operator has stated the answer, and
       honouring it costs no taxonomy lookup at all.
    2. **Direct membership**: the queried taxid is itself in a family's
       `taxonomic_scope`. Short-circuits before any lineage lookup (no
       network/subprocess call), so this stays free for the common case.
    3. **Lineage membership**: the taxid's NCBI Taxonomy ancestor lineage
       contains a taxid in a family's `taxonomic_scope`. Most families declare
       a broad scope (a subphylum/subclass taxid) expecting it to cover every
       descendant species; lineage matching is what makes that work. Before it
       existed, an audit found 50 of 61 curated records fell through to rule 5.
    4. **Phylum fallback**: nothing above matched, but the query taxon's own
       phylum is known and the database has families for it. Only that phylum's
       families are searched. A MAT locus family curated in Basidiomycota
       cannot be the answer for an Ascomycota genome, so searching it is pure
       cost: it makes the run slower AND adds cross-phylum candidates for the
       polish stage to grind through. A live case -- taxid 294748, lineage
       [131567, 2759, 33154, 4751, 451864, 4890, 716545, 147537, 3239874,
       2916678, 766764, 5475, 5476] -- intersects NO curated family's scope, so
       before this rule it routed to all 19 families across all three phyla and
       one small yeast genome ran past 24 minutes without finishing.
    5. **Exhaustive**: every family. Reached when there is no taxid at all, or
       the phylum cannot be determined, or the determined phylum has no curated
       families (returning that phylum's EMPTY family set instead would silently
       detect nothing, which is worse than searching too much).

    **How the phylum is determined.** Not from a hardcoded phylum-name-to-taxid
    table, and not by testing which phylum's families have a scope taxid in the
    lineage -- the latter is exactly the test rule 3 just failed, so by
    construction it can never succeed here. Instead `phylum_name_resolver`
    returns the NAME of the lineage's phylum-rank ancestor, which is compared
    against `FamilyKey.phylum`; both sides are then NCBI scientific names, and
    `available_phyla` guarantees the database side is read from disk. The
    default resolver (`db.taxonomy.default_lineage_phylum_name`) parses the very
    same cached `efetch db=taxonomy` document the default lineage resolver
    parses, at the same URL, so this fallback adds NO network call to a run that
    already fetched the lineage.

    Every resolver failure (network error, unknown taxid) degrades to the next
    rule rather than raising, so a taxonomy outage makes detection slower, never
    broken.
    """
    if phylum is not None:
        return RoutingDecision(
            families=[f for f in families if f.key.phylum == phylum],
            routing_mode="explicit_phylum",
            phylum=phylum,
        )
    if taxid is None:
        return RoutingDecision(families=list(families), routing_mode="exhaustive")

    direct = [f for f in families if taxid in f.taxonomic_scope]
    if direct:
        return RoutingDecision(families=direct, routing_mode="direct")

    try:
        ancestors = set(lineage_taxids_resolver(taxid))
    except Exception:
        ancestors = set()
    lineage_matched = [f for f in families if ancestors & set(f.taxonomic_scope)]
    if lineage_matched:
        return RoutingDecision(families=lineage_matched, routing_mode="lineage")

    try:
        query_phylum = phylum_name_resolver(taxid)
    except Exception:
        query_phylum = None
    if query_phylum:
        in_phylum = [f for f in families if f.key.phylum == query_phylum]
        if in_phylum:
            return RoutingDecision(
                families=in_phylum, routing_mode="phylum_fallback", phylum=query_phylum
            )

    return RoutingDecision(families=list(families), routing_mode="exhaustive")
