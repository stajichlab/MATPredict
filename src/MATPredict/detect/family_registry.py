"""Load order.yml families and route them to a taxid via taxonomic_scope."""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
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


DEFAULT_MAX_HOMOTHALLIC_SEPARATION_BP = 20_000
"""How close two opposite-idiomorph genes must be to be called one locus.

Both idiomorphs present in one genome is the homothallic architecture, not an
error: Syzygites megalocarpus encodes both HMG transcription factors, each
flanked by its own intact gene with the other flank pseudogenised (Idnurm 2011,
doi:10.1128/EC.05149-11; Schulz et al. 2016). Calling that a mistake would make
homothallism undetectable.

But two idiomorph genes merely landing in one cluster is not enough. The
clustering gap for Mucoromycota is 50 kb, so an unrelated pair can be grouped:
measured on Syzygites sp. MES_3091, a real annotated sexM and a 285 bp tblastn
sexP fragment 47 kb apart were grouped as one "locus". A real homothallic pair
spans genes, not tens of kb of nothing.

PROVISIONAL, curator-set 2026-09-20 as a starting value. Revise it from
observed separations in confirmed homothallic loci, not by intuition.
"""

DEFAULT_MAX_PLAUSIBLE_LOCUS_SPAN_BP = 200_000
"""The widest a locus of this family is expected to be, in bp.

A FLAG, never a filter. Curator ruling, J. Stajich, 2026-09-20: large MAT
loci are real and must still be found -- consistent with unpublished findings
by a former graduate student -- so a wide call is reported and marked, not
dropped.

Measured basis, 283-genome BFD Mucoromycota sweep (687 loci): median span
13,565 bp, p90 57,085, p99 105,222, max 154,355. Ground-truth loci in
`testset/Zygo` span 6,795-13,089 bp. At 200 kb this flag fires on NOTHING
already observed; it is a guard-rail against a runaway cluster rather than a
filter on present output.

Why not the 120 kb first considered: it would have flagged 3 of 687, one of
them a Blakeslea trispora call with six genes at 87.2% identity sitting
473 bp over the line. Blakeslea is itself a curated reference, so that
identity is partly same-species, but a six-gene call is not what a bound of
this kind should be catching.

PROVISIONAL. Revise from observed spans in confirmed loci, not by intuition.
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
    #: Curated gene name -> the roster's canonical name for it. Built from each
    #: gene entry's optional `aliases:` list, plus every canonical name mapping
    #: to itself, so a caller can look any name up without a special case.
    #:
    #: Exists because a curated record deposits the gene name its PUBLICATION
    #: used, and several of those can be one gene. In 40410_jec20_MAT_a,
    #: MFa1/MFa2/MFa3 are three byte-identical 42 aa proteins; no ranking rule
    #: can separate them, so they are one roster gene with three aliases. The
    #: records keep their deposited names; only the roster collapses.
    gene_aliases: dict[str, str] = field(default_factory=dict)
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
    max_homothallic_separation_bp: int = DEFAULT_MAX_HOMOTHALLIC_SEPARATION_BP
    """How close opposite-idiomorph genes must be to call one homothallic locus.

    Curation data like `max_cluster_gap_bp`: how compact a homothallic MAT
    region is, is a property of the clade's locus architecture. See
    `DEFAULT_MAX_HOMOTHALLIC_SEPARATION_BP`.
    """
    max_plausible_locus_span_bp: int = DEFAULT_MAX_PLAUSIBLE_LOCUS_SPAN_BP
    """The widest a locus of this family is expected to be, in bp.

    Curation data like `max_cluster_gap_bp`, and a FLAG rather than a filter:
    a call wider than this is still reported, marked so a reader can see it.
    See `DEFAULT_MAX_PLAUSIBLE_LOCUS_SPAN_BP` for the measured basis.
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
    #: Why a resolver failed, when one did. Degrading is by design, but doing
    #: it invisibly let one taxid route three different ways in one pilot
    #: with no report saying why. None when every lookup succeeded.
    routing_error: str | None = None


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


def _gene_alias_map(genes: list) -> dict[str, str]:
    """{curated gene name -> canonical roster name}, including each canonical
    name mapping to itself so callers need no special case.

    Tolerates a gene entry that is not a mapping (some order.yml fixtures
    declare a bare name); such an entry contributes nothing.
    """
    mapping: dict[str, str] = {}
    for gene in genes or []:
        if not isinstance(gene, dict) or "name" not in gene:
            continue
        canonical = gene["name"]
        mapping[canonical] = canonical
        for alias in gene.get("aliases") or []:
            mapping[alias] = canonical
    return mapping


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
                    # Built eagerly, so it must tolerate whatever `genes` holds:
                    # this loader passes the list through unnormalised and some
                    # callers declare a gene as a bare name rather than a
                    # mapping. A non-mapping entry simply has no aliases.
                    gene_aliases=_gene_alias_map(locus.get("genes", [])),
                    min_idiomorph_margin=locus.get(
                        "min_idiomorph_margin", DEFAULT_MIN_IDIOMORPH_MARGIN
                    ),
                    max_homothallic_separation_bp=locus.get(
                        "max_homothallic_separation_bp",
                        DEFAULT_MAX_HOMOTHALLIC_SEPARATION_BP,
                    ),
                    max_plausible_locus_span_bp=locus.get(
                        "max_plausible_locus_span_bp",
                        DEFAULT_MAX_PLAUSIBLE_LOCUS_SPAN_BP,
                    ),
                )
            )
    return families


def derive_max_cluster_gap(
    families: list[Family], routing_mode: str | None = None
) -> int:
    """The clustering gap for a run over `families`: the MAXIMUM of their gaps,
    EXCEPT on a failed route, where the default stands instead.

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

    `routing_mode` qualifies all of the above, added 2026-09-21. The "maximum
    is the safe end" argument holds only among families that could plausibly
    BE this genome's locus -- which is exactly what a matched route
    establishes and a failed one does not. On `phylum_fallback` and
    `exhaustive` the run searches every family in the phylum precisely
    BECAUSE nothing matched, so inheriting the widest outlier's gap is not a
    conservative choice, it is an arbitrary one.

    This became load-bearing when the curator set the Tremellales `MAT` gap
    to 120 kb (the Cryptococcus MAT locus really does span ~104 kb with a
    74.9 kb internal gene gap, measured from AF542531.2/AF542530.2). Without
    this qualifier that one locus would cluster EVERY unrouted Basidiomycota
    genome at 120 kb -- roughly 1,299 BFD genomes even after order-level
    scoping, since Boletales, Polyporales, Sporidiobolales, Trichosporonales,
    Pucciniales and Cantharellales have no curated record between them.

    `explicit_phylum` is deliberately NOT capped: `--phylum` is an operator
    assertion about the query, not a failed lookup, and the validated
    Mucoromycota workflow runs `--phylum Mucoromycota` and depends on that
    locus's curated 50 kb gap. Capping it would silently halve that and break
    the 23/23 ground-truth result.

    Omitting `routing_mode` keeps the old behaviour exactly, so callers that
    do not know the mode are unaffected.
    """
    if routing_mode in ("phylum_fallback", "exhaustive"):
        return DEFAULT_MAX_CLUSTER_GAP_BP
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


#: `load_record_families` results, keyed by db_root, with the directory stamp
#: they were built from. `run_pipeline` calls that function once per genome and
#: it parses every curated metadata.yaml: 465 ms against the live 101-record
#: database, versus 2.5 ms to stat the same files. Over a 3,174-genome rollout
#: the repeated parse is roughly 27 minutes. Invalidated by content hash.
_RECORD_FAMILIES_CACHE: dict[Path, tuple[str, dict[str, "FamilyKey"]]] = {}


def clear_record_families_cache() -> None:
    """Drop the cached record->family indexes. For tests, and for any caller
    that edits `db/` in-process and wants the next read to be honest."""
    _RECORD_FAMILIES_CACHE.clear()


def _db_stamp(db_root: Path) -> str:
    """A content hash of every curated metadata.yaml under `db_root`.

    A real hash, not a (size, mtime) heuristic. The heuristic was tried first
    and is unsafe: `st_mtime_ns` granularity is filesystem-dependent, and on
    this cluster's `/scratch` five rapid rewrites of one file produced only
    TWO distinct mtimes. A curator changing a `locus_name` from `HD` to `PR`
    -- same byte length, so same file size and same file count -- would then
    be served a stale index, and every hit in the run would be attributed to
    the wrong family. The repo's own test suite caught exactly that case.

    Hashing is cheap enough that the heuristic bought nothing: 5.4 ms to hash
    all 101 curated records against 466 ms to parse them, an 86x margin. The
    path is hashed alongside the bytes so a rename, or a record moved between
    phyla, invalidates too.
    """
    digest = hashlib.blake2b(digest_size=16)
    for path in sorted(db_root.glob("*/*/*/metadata.yaml")):
        digest.update(str(path.relative_to(db_root)).encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


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
    stamp = _db_stamp(db_root)
    cached = _RECORD_FAMILIES_CACHE.get(db_root)
    if cached is not None and cached[0] == stamp:
        return dict(cached[1])

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
    _RECORD_FAMILIES_CACHE[db_root] = (stamp, dict(index))
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

    # Exact-taxid and lineage matches are UNIONED, not cascaded. Returning at
    # the first tier that matched anything let a SPECIES-scoped family shadow a
    # GENUS-scoped one for the same organism, and silently drop the latter.
    #
    # Measured on Schizosaccharomyces pombe (4896), which is how this was
    # found: `mat2` and `mat3` are scoped [4896] and matched exactly, so the
    # function returned those two and never reached the lineage tier -- where
    # `mat1`, scoped [4895] (the genus Schizosaccharomyces, 4896's parent),
    # would have matched. `mat1` is the ACTIVE mating-type locus; mat2 and
    # mat3 are the silent cassettes. Detection therefore looked like it worked
    # on S. pombe while never searching for the locus that actually determines
    # mating type. The leave-one-out benchmark caught it as an unexplained miss
    # with a reference still present.
    #
    # A genus-level scope is not weaker evidence than a species-level one; it
    # is a curator's statement about a different breadth. Nothing justifies one
    # suppressing the other. `routing_mode` still reports `direct` when an
    # exact match contributed, so the stronger signal stays visible.
    errors: list[str] = []
    direct = [f for f in families if taxid in f.taxonomic_scope]
    try:
        ancestors = set(lineage_taxids_resolver(taxid))
    except Exception as exc:
        ancestors = set()
        errors.append(f"lineage lookup failed: {type(exc).__name__}: {str(exc)[:160]}")
    lineage_matched = [f for f in families if ancestors & set(f.taxonomic_scope)]
    if direct or lineage_matched:
        seen: set = set()
        merged = []
        for f in direct + lineage_matched:
            if f.key not in seen:
                seen.add(f.key)
                merged.append(f)
        return RoutingDecision(
            families=merged, routing_mode="direct" if direct else "lineage",
            routing_error="; ".join(errors) or None,
        )

    try:
        query_phylum = phylum_name_resolver(taxid)
    except Exception as exc:
        query_phylum = None
        errors.append(f"phylum lookup failed: {type(exc).__name__}: {str(exc)[:160]}")
    if query_phylum:
        in_phylum = [f for f in families if f.key.phylum == query_phylum]
        if in_phylum:
            return RoutingDecision(
                families=in_phylum, routing_mode="phylum_fallback", phylum=query_phylum,
                routing_error="; ".join(errors) or None,
            )

    return RoutingDecision(families=list(families), routing_mode="exhaustive",
                           routing_error="; ".join(errors) or None)
