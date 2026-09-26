"""Orchestrates routing -> search -> clustering -> polishing -> scoring -> tiering -> idiomorph assignment.

Search is a two-stage *localize-then-polish* flow (see
docs/superpowers/specs/2026-09-17-mat-detection-search-localization-design.md):

* **Genome-only path** (no predicted proteome): one batched genome-wide
  `tblastn` localization (`search_localize`) over every routed family's
  reference proteins, clustered once by `cluster_hits`. Every gene that
  `tblastn` localized in a cluster is then polished.
* **Fast-path windowed rescue** (a proteome was supplied): `search_fast_path`
  runs as before; a cluster/family that has a foothold but is still missing one
  of its own `core_MAT` genes polishes that one gene directly against a window
  padded around the *existing* cluster's span -- a rough location already
  exists.
* **Fast-path localization rescue** (a proteome was supplied): ONE batched
  `search_localize` (`tblastn`) call, never one per family, covering both
  * every routed family with NO fast-path hit at all (no foothold, so no window
    to polish against), for its whole gene set, and
  * every family with a PARTIAL foothold, restricted to the specific `core_MAT`
    genes the proteome did not place in some *individual* fast-path cluster of
    that family (see `_localization_rescue_targets`). Eligibility is keyed per
    `(family, gene, cluster)`, never per `(family, gene)`: a family can have
    more than one real, independent locus in one genome (tetrapolar species
    with unlinked loci; homothallic/heterothallic switching-cassette species),
    so a gene found in cluster A is no evidence that cluster B has it.
  This is the blind spot the whole pipeline exists to close: a short
  pheromone-precursor gene that a supplied genome annotation simply does not
  contain cannot be found by searching that annotation, and it need not sit
  inside the narrow window around whatever else of its family WAS annotated.
  The rescued `tblastn` hits are clustered together with the fast-path hits and
  are polished exactly like any other localized cluster. Both rescues run: the
  windowed one refines the neighbourhood already known, the localization one
  looks genome-wide.

Polish eligibility is therefore decided per (cluster, family, gene), never by a
single global "am I in genome-only mode" flag: a gene is polished when it was
localized by `tblastn` (whichever path produced that localization) or when it is
one of its family's own `core_MAT` genes still missing from that cluster -- but
never when that same cluster already holds a non-localized (annotated) hit for
that gene, which contains the guard/chaining mismatch described on
`_RescueScope.accepts` and in the polish loop below.

Polishing runs BOTH `exonerate --refine region` and `miniprot` against the same
padded window and classifies the pair (`polish.classify`) into one of
`polished_agree` / `polished_disagree` / `polished_single` / `unpolished`.
Only `unpolished` (neither tool produced a model, but the raw localization hit
stands) affects confidence tiering, capping the family at Medium -- exactly the
effect the retired "relaxed exonerate second pass" used to have. Tool agreement
itself is reported but never consulted by `tiering.assign_tier`.

Fragmented assemblies (spec section 6) are handled, OPT-IN only, to the extent
described in `_fragmented_family_segments`: a family whose expected `core_MAT`
genes are split across clusters on *different contigs*, with no single cluster
carrying them all, can be reported as one multi-segment call with
`fragmented=True` (which downgrades its confidence tier by one, per
`tiering.assign_tier`) -- but only when the caller passes
`allow_cross_contig_fragments=True`. `contig_edge_distance` is populated per
segment when the genome FASTA can be read; it is left as None otherwise. What
is deliberately NOT attempted here is reconstructing locus order or the
intervening sequence across segments -- a multi-segment call reports the
segments it found, nothing more.

The DEFAULT is `allow_cross_contig_fragments=False`: cross-contig merging is
OFF. Curator's ruling, 2026-09-20, after Basidiomycota order testing found a
real Cryptococcus deneoformans MAT locus (both genes at 100% identity)
reported as one merged call whose top-level contig/start/end named a
chromosome holding almost none of its own evidence, giving
`idiomorph:undetermined` and `confidence: low` for a strain whose mating type
was not in doubt. With the merge disabled, each contig's own partial evidence
is still reported -- by the ordinary per-cluster loop, not this function --
as separate, honest, correctly-coordinated calls.
"""
from __future__ import annotations

import json
import logging
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from MATPredict.detect.clustering import GeneCluster, cluster_hits
from MATPredict.detect.family_registry import (
    Family,
    FamilyKey,
    RoutingDecision,
    derive_max_cluster_gap,
    load_all_families,
    load_record_families,
    route,
)
from MATPredict.detect.idiomorph import (
    LOCUS_CLASS_MAT,
    LOCUS_CLASS_PARTIAL,
    apply_partial_locus,
    idiomorph_candidates,
    idiomorph_margin_from_vote,
    is_partial_strength,
    DEFAULT_MIN_OVERLAP_FRACTION,
    IdiomorphResolution,
    assign_idiomorph,
    classify_locus,
    evidenced_idiomorphs,
    resolve_idiomorph_overlaps,
)
from MATPredict.detect.polish import (
    STATUS_AGREE,
    STATUS_DISAGREE,
    STATUS_SINGLE,
    STATUS_NOT_POLISH_CANDIDATE,
    STATUS_UNPOLISHED,
    PolishModel,
    PolishOutcome,
    classify,
)
from MATPredict.db.taxonomy import default_genetic_code
from MATPredict.detect.reference_fasta import searchable_genes_by_family
from MATPredict.detect.scoring import FamilyScore, is_ambiguous, score_cluster
from MATPredict.detect.search import (
    SearchHit,
    drop_low_quality_hits,
    polish_with_exonerate,
    polish_with_miniprot,
    search_fast_path,
    search_localize,
)
from MATPredict.detect.tiering import assign_tier, cap_at_medium

logger = logging.getLogger(__name__)

#: out_path strings already warned about by `_write_evidence_diagnostics` in
#: this process, so a bad path (typo, missing directory) logs one warning per
#: run instead of spamming one per cluster/family row.
_DIAGNOSTICS_WRITE_FAILURES_LOGGED: set[str] = set()


@dataclass(frozen=True)
class GeneEvidence:
    """Per-gene detail behind a detection, kept so the report/GFF3 writers can
    emit identity, coverage, coordinates, role and the curated record matched
    (spec section 8) rather than just a gene name."""

    gene_name: str
    role: str
    contig: str
    start: int
    end: int
    strand: str
    identity: float
    coverage: float | None
    reference_record_id: str
    method: str
    status: str = "polished_agree"  # one of polish.STATUS_* -- default only
    # for backward-compatible construction in existing tests that predate
    # this field; real pipeline code below always passes a real, computed
    # status, never this default.
    alternate_model: dict | None = None  # the OTHER tool's model
    # (contig/start/end/strand/exons/identity/method), populated only when
    # status == polish.STATUS_DISAGREE, else None.
    exons: tuple[tuple[int, int], ...] | None = None
    # The canonical polished model's own real exon structure (1-based,
    # fully-closed genomic spans), populated only when this evidence came
    # from a polished `PolishModel` with a non-empty `exons` list. **Order
    # is ASCENDING GENOMIC COORDINATE, regardless of strand** -- this
    # matches how `search.py`'s `polish_with_exonerate`/`polish_with_miniprot`
    # actually build `PolishModel.exons` (`sorted(..., key=lambda e:
    # int(e[3]))`, i.e. by genomic start, for both strands), which is NOT
    # the same convention as `db/ncbi_client.py`'s `CdsStructure` (already
    # in transcript order, descending for minus strand). A caller that
    # needs a spliced transcript (e.g.
    # `detect.benchmark._extract_translated_gene`) must reverse this list
    # itself for a minus-strand gene before concatenating. A raw, unpolished
    # `SearchHit` genuinely has no exon structure to offer (its
    # `identity`/`coverage` come from one ungapped or splice-naive
    # alignment), so this stays `None` for that fallback path -- callers
    # that need a protein sequence for a `None`-exons gene fall back to
    # naive single-span translation for that gene, which is inherently
    # approximate for a real multi-exon gene reported this way.


@dataclass(frozen=True)
class LocusSegment:
    contig: str
    start: int
    end: int
    contig_edge_distance: int | None = None


@dataclass(frozen=True)
class DetectionResult:
    family_key: FamilyKey
    contig: str
    start: int
    end: int
    confidence: str
    idiomorph: str
    ambiguous_with: list[FamilyKey]
    genes_found: list[str]
    genes_missing: list[str]
    fragmented: bool
    genes_not_searchable: list[str] = field(default_factory=list)
    segments: list[LocusSegment] = field(default_factory=list)
    gene_evidence: list[GeneEvidence] = field(default_factory=list)
    reference_records: list[str] = field(default_factory=list)
    locus_class: str = "mat_locus"
    """WHAT was found: `mat_locus`, `homothallic_candidate`,
    `idiomorph_gene_only` or `flanking_gene_only`. Orthogonal to
    `detection_pass`, which says HOW it was admitted. See
    `idiomorph.classify_locus`; `idiomorph_gene_only` and `flanking_gene_only`
    are kept deliberately as leads, not discarded."""
    span_exceeds_plausible_bound: bool = False
    """Is this call wider than its family's `max_plausible_locus_span_bp`?

    A flag, never a filter -- the curator ruled on 2026-09-20 that wide loci
    are real and must still be reported. Present on EVERY result, so a reader
    can tell "not flagged" from "this build predates the field". See
    `span_exceeds_plausible_bound` (the function) and
    `family_registry.DEFAULT_MAX_PLAUSIBLE_LOCUS_SPAN_BP`.
    """
    detection_pass: str = "strict"
    """Which pass admitted this locus: `strict` or `relaxed`.

    Present on EVERY result, not only relaxed ones, so a reader can tell
    "this was a strict call" from "this build predates the field". A relaxed
    call cleared only the evidence floor (>=2 distinct genes, >=1 core_MAT)
    and not the fraction floor, so its tier is capped and anyone who wants
    strict-only output can filter on this field.
    """
    #: Every idiomorph this cluster has evidence for, best score first. Emitted
    #: alongside the scalar `idiomorph` so an exact tie, or a narrow win, is
    #: visible rather than collapsed. Curator's ruling 2026-09-21: on a tie
    #: "report both instead of worrying about getting it right".
    idiomorph_candidates: list[dict] = field(default_factory=list)
    idiomorph_margin: float | None = None
    """Identity points separating the winning idiomorph gene from the loser.

    `None` when no idiomorph resolution was needed for this locus, which is
    the ordinary case for a family whose genes are not mutually exclusive.
    When several resolutions contributed, this is the NARROWEST of them --
    the call is only as trustworthy as its weakest step. Below the locus's
    `min_idiomorph_margin` the confidence tier is capped; the number is
    reported either way so a close call is never mistaken for a clean one.
    """
    polished_genes: int = 0
    """How many of this locus's genes produced a real polished gene model.

    A gene counted here was modelled by `exonerate --refine` or `miniprot`
    (`polished_agree`, `polished_disagree` or `polished_single`). A raw tblastn
    HSP that no tool could turn into a gene (`unpolished`) does NOT count, and
    neither does a gene evidenced directly from an annotation
    (`not_polish_candidate`) -- the question this answers is "did a gene model
    survive here", and an annotated gene was never asked.

    This is the sharpest single discriminator measured on the 2026-09-22
    Pezizomycotina panels (46,647 lineage-routed loci):

        polished genes    loci        high   mean identity
                     0  39,838 (85%)     0           32.8%
                     1   1,161 ( 2%)     0             --
                    2+   5,648 (12%)  2,957           61.9%

    EVERY high-confidence call has two or more; not one of the 41,000 loci
    below that bar is high-confidence, and their mean identity sits in the
    twilight zone where alignment stops implying homology. They are not a
    weaker tail of the real signal, they are scattered HMG-box and alpha-box
    paralogs -- 150-350 bp HSPs at 8-28% query coverage.
    """
    idiomorph_resolutions: list[IdiomorphResolution] = field(default_factory=list)
    """Every overlapping idiomorph pair collapsed for this locus.

    Kept in full because a novel or hybrid locus is exactly what a low margin
    might indicate, and the curator asked that the ambiguity be reported
    rather than silently resolved away. Also the observations a future
    recalibration of the overlap threshold needs.
    """


@dataclass(frozen=True)
class NotDetectedFamily:
    """An attempted family that produced no reported locus, and why.

    The spec requires a sub-floor outcome to be reported as "not detected",
    explicitly listing which families were attempted and why each fell short,
    never silently omitted.
    """

    family_key: FamilyKey
    reason: str
    best_fraction_found: float
    genes_found: list[str] = field(default_factory=list)
    genes_missing: list[str] = field(default_factory=list)
    genes_not_searchable: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class DetectionOutcome:
    """Everything one `matpredict detect` run concluded: the loci it called,
    the attempted families it could not call, and the full attempted set."""

    results: list[DetectionResult]
    not_detected: list[NotDetectedFamily] = field(default_factory=list)
    families_attempted: list[FamilyKey] = field(default_factory=list)
    suppressed_unpolished: int = 0
    """How many built loci were withheld for having fewer than
    MIN_POLISHED_GENES polished genes.

    Reported as a number so the suppression is visible rather than silent: a
    genome that drops from 7 loci to 1 says so. The per-candidate detail is
    already on disk in the evidence-diagnostics stream, so nothing withheld
    here is lost.
    """
    suppressed_loci: list[DetectionResult] = field(default_factory=list)
    """The withheld loci themselves, so a bar loss can be audited against a
    known locus. Without their coordinates a withheld call at the right place
    looks exactly like no call at all."""
    # Which `family_registry.route` rule chose `families_attempted` (see
    # `RoutingDecision`). Carried all the way into the detection report because
    # it is what tells a reader whether this run's not-detected entries are
    # evidence of absence or just the exhaustive fallback searching phyla the
    # query could not belong to. None when the outcome was built without
    # routing information (a direct `run_pipeline` call in a test).
    routing_mode: str | None = None
    #: A taxonomy lookup that failed and so widened the routing; None when
    #: every lookup succeeded. See `RoutingDecision.routing_error`.
    routing_error: str | None = None
    #: The translation table this run used, and why it fell back to table 1
    #: when the lookup failed (None when it did not fail).
    genetic_code: int | None = None
    genetic_code_error: str | None = None


def _missing_core_genes(cluster: GeneCluster, family: Family) -> set[str]:
    """A single family's own core_MAT genes not yet found among that same
    family's own hits in this cluster -- never checked against another
    family's hits or expected genes (see Task 6's cross-family pooling bug)."""
    core_genes = {g["name"] for g in family.genes if g["role"] == "core_MAT"}
    found_genes = {h.gene_name for h in cluster.hits if h.family_key == family.key}
    return core_genes - found_genes


#: A polished model must cover at least this fraction of the reference protein
#: it was built from to count toward a `homothallic_candidate` of two unrelated
#: genes. Measured 2026-09-25: at 0.5 the rule fires on 0 of 113 would-be loci
#: in heterothallic-rich panels and keeps both Hydnotrya loci; at 0.6 it also
#: loses those. PROVISIONAL, like the other curator-tunable thresholds.
HOMOTHALLIC_MIN_MODEL_COVERAGE = 0.5


def full_length_models(
    polish_by: dict,
    cluster_ids: set[int],
    family_key: FamilyKey,
    record_lengths: dict[tuple[str, str], int],
    cross_matched: set[str],
    min_coverage: float = HOMOTHALLIC_MIN_MODEL_COVERAGE,
) -> frozenset[str]:
    """Genes of this family, in these clusters, that are full-length models.

    A gene qualifies when a polishing tool modelled it (agree, disagree or
    single), the canonical model covers >= `min_coverage` of the reference
    protein it was built from (`record_lengths[(record_id, gene)]`, in aa),
    and the gene took no part in an idiomorph cross-match resolution here
    (`cross_matched`: every winner and loser). Used only by the relaxed
    homothallic rule: a truncated remnant, or an HMG gene that beat another
    HMG gene at one position, is not evidence of a second idiomorph.
    """
    out = set()
    for (cluster_id, key, gene_name), outcome in polish_by.items():
        if key != family_key or cluster_id not in cluster_ids or gene_name in cross_matched:
            continue
        if outcome.status not in (STATUS_AGREE, STATUS_DISAGREE, STATUS_SINGLE):
            continue
        model = outcome.canonical
        reference = record_lengths.get((model.reference_record_id, gene_name))
        if not reference:
            continue
        aa = sum(e.end - e.start + 1 for e in model.exons) // 3
        if aa >= min_coverage * reference:
            out.add(gene_name)
    return frozenset(out)


def _curated_protein_lengths(
    db_root: Path,
    families: list[Family],
    record_families: dict[str, FamilyKey],
    exclude_record_ids: frozenset[str] = frozenset(),
    by_record: dict[tuple[str, str], int] | None = None,
) -> dict[tuple[FamilyKey, str], int]:
    """`(family_key, gene_name)` -> the LONGEST curated reference protein, in aa.

    When `by_record` is given it is also filled with `(record_id, gene_name)`
    -> that record's protein length, for `full_length_models`.

    Keyed per `(phylum, locus_name)` family, never by bare gene name: gene names
    are reused across families (`sla2` is both `Ascomycota:MATsc`'s and
    `Ascomycota:MATyl`'s), so one family's entry must never answer for another's.
    See `_short_orf_genes` for why the MAXIMUM is the right summary.

    Reads db/**/proteins.faa (matching gff_export.write_proteins_fasta's header
    form `>{record_id}|gene_index={n}|name={name}|role={role}`, the same raw
    per-record files reference_fasta.py concatenates -- not its rewritten output).
    """
    expected_by_family = {f.key: {g["name"] for g in f.genes} for f in families}
    # Curated records carry the gene name their publication deposited; the
    # roster may have collapsed several of those onto one canonical name (see
    # `Family.gene_aliases`). Resolve through it so a collapsed gene's length
    # is found rather than silently skipped by the `expected_by_family` test.
    aliases_by_family = {f.key: f.gene_aliases for f in families}
    longest: dict[tuple[FamilyKey, str], int] = {}

    for faa in db_root.glob("*/*/*/proteins.faa"):
        if faa.relative_to(db_root).parts[0] == "candidates":
            continue
        for chunk in faa.read_text().split(">")[1:]:
            header, _, seq = chunk.partition("\n")
            head_fields = header.split("|")
            record_id = head_fields[0]
            # A withheld record must not reach ANY db-derived quantity, not
            # just the query set. Its protein length sets the polish window's
            # padding and its length decides the short-ORF/unsearchable split,
            # so leaving it in would let a held-out answer shape the run that
            # is supposed to be blind to it.
            if record_id in exclude_record_ids:
                continue
            family_key = record_families.get(record_id)
            if family_key is None or family_key not in expected_by_family:
                continue
            parts = dict(p.split("=", 1) for p in head_fields[1:] if "=" in p)
            name = aliases_by_family.get(family_key, {}).get(
                parts.get("name"), parts.get("name")
            )
            if name not in expected_by_family[family_key]:
                continue
            length = len(seq.strip().replace("\n", ""))
            key = (family_key, name)
            longest[key] = max(longest.get(key, 0), length)
            if by_record is not None:
                by_record[(record_id, name)] = max(by_record.get((record_id, name), 0), length)
    return longest


def _short_orf_genes(
    db_root: Path,
    families: list[Family],
    record_families: dict[str, FamilyKey],
    floor_aa: int,
    exclude_record_ids: frozenset[str] = frozenset(),
) -> dict[FamilyKey, set[str]]:
    """Per-family gene names whose BEST curated reference protein is shorter than floor_aa.

    Scoped per `(phylum, locus_name)` family and keyed on the MAXIMUM curated
    length seen for that gene *within that family*, for two reasons:

    * Global scoping is wrong because gene names are reused across families
      (`sla2` is both `Ascomycota:MATsc`'s flanking gene and `Ascomycota:MATyl`'s),
      so one family's short entry must not condemn another family's gene.
    * The minimum (or first-seen) length is wrong because the curated database
      legitimately contains fragments: `Ascomycota:MATsc` holds a 22-aa `cha1`
      fragment alongside normal-length curated proteins, and taking the shortest
      made the pipeline report the perfectly searchable `cha1` gene as "not
      searchable by this method" -- the exact false claim this feature exists to
      prevent. Maximum is chosen over median because the question being asked is
      "is there ANY curated protein long enough to search with?", and a single
      full-length example is enough to make the gene searchable.

    The per-family/per-gene length scan itself lives in
    `_curated_protein_lengths`, shared with the polish-window padding helper.
    """
    longest = _curated_protein_lengths(
        db_root, families, record_families, exclude_record_ids)

    short_by_family: dict[FamilyKey, set[str]] = {}
    for (family_key, name), length in longest.items():
        if length < floor_aa:
            short_by_family.setdefault(family_key, set()).add(name)
    return short_by_family


@dataclass(frozen=True)
class EvidenceFloor:
    """Minimum evidence a family must clear in a cluster before it is admitted
    to the (expensive, per-gene, two-subprocess-per-tool) Stage 2 polish loop.

    The defaults are the CURATOR'S RULING: a family needs **>=2 distinct genes
    in the cluster, at least one of them with role `core_MAT`**. A flanking
    partner is deliberately NOT required. Requiring one would systematically
    miss fragmented assemblies, where the MAT locus is real but its flank sits
    on another contig -- exactly the genomes this tool is most needed for. In
    the curator's words: "polish trigger could be another gene not necessarily
    flank - at least at start until we have a bigger collection of examples of
    fragments."

    `min_identity` stays `None` on purpose. The ruling was about gene count and
    role, not identity, and no identity cutoff has been measured against real
    divergent MAT proteins. An unmeasured identity floor would silently drop
    true positives.

    These defaults previously reproduced the retired `_families_with_a_foothold`
    behavior (any single hit, any role, any identity), with tuning deferred
    until routing was narrowed. That deferral is RESOLVED: phylum-aware routing
    now sends a query to its own phylum's families instead of the exhaustive
    fallback that produced the noisy calibration example, so the floor no longer
    has to absorb cross-phylum noise on its own.

    Note that a default of `min_hits=1, require_core_role=False` is still what
    the evidence diagnostics enumerate with (`_DIAGNOSTICS_CANDIDATE_FLOOR`):
    the diagnostics record every candidate, admitted or not, and remain the
    dataset for any future tightening (e.g. of `min_identity`).
    """

    #: Minimum number of DISTINCT genes (by `SearchHit.gene_name`), never raw
    #: hit/HSP objects. `search.py`'s tblastn loop appends one `SearchHit` per
    #: HSP line with no per-gene deduplication, so one real gene matched by N
    #: curated reference records (routine as a locus's own curated database
    #: grows) produces N hits -- counting raw hits here would make this floor's
    #: effective strictness silently drift with curation density rather than
    #: with actual evidence. See `_families_meeting_evidence_floor`.
    min_hits: int = 2
    min_identity: float | None = None
    require_core_role: bool = True


#: How many polished gene models a cluster must carry to be REPORTED as a locus.
#: Curator's ruling, 2026-09-22, having been shown that loci with exactly one
#: polished gene produce zero high-confidence calls: "yes require more than 1".
#: Measured cost of the bar on the 2026-09-22 Pezizomycotina panels: 46,647
#: reported loci fall to 5,648 (12.1%) and ALL 2,957 high-confidence calls
#: survive. At >=1 the figure would be 6,809 and the same 2,957 -- so the
#: stricter bar removes a further 1,161 loci for nothing.
MIN_POLISHED_GENES = 2

#: The permissive floor the evidence diagnostics enumerate candidates with:
#: every family with >=1 own hit in the cluster, of any gene, any role, any
#: identity. This is the retired `_families_with_a_foothold` behavior, and it
#: is spelled out here rather than as `EvidenceFloor()` because `EvidenceFloor`'s
#: own defaults are now the strict curator-ruled floor -- writing `EvidenceFloor()`
#: for the diagnostics would shrink the calibration dataset to just the rows
#: that were admitted anyway, which is precisely the information it must not lose.
_DIAGNOSTICS_CANDIDATE_FLOOR = EvidenceFloor(
    min_hits=1, min_identity=None, require_core_role=False
)


def _own_live_hits(cluster, family_key):
    """This family's non-superseded hits in `cluster`.

    Used for the idiomorph vote, which ranks evidence rather than testing
    presence. A superseded hit is the losing half of a resolved cross-match and
    must not vote for its own idiomorph.
    """
    return [
        h for h in cluster.hits
        if h.family_key == family_key and h.superseded_by is None
    ]


def _families_meeting_evidence_floor(
    cluster: GeneCluster, families: list[Family], floor: EvidenceFloor
) -> list[Family]:
    """Families whose OWN hits in this cluster clear `floor` -- generalizes
    the retired `_families_with_a_foothold` (which was exactly
    `_families_meeting_evidence_floor(cluster, families,
    _DIAGNOSTICS_CANDIDATE_FLOOR)`)."""
    admitted = []
    for family in families:
        # Superseded hits are excluded before anything is counted. Such a hit
        # and its winner hit the SAME locus gene under two mutually exclusive
        # idiomorph names, so counting it would let a cluster whose only
        # evidence is one HMG gene clear a bar that asks for two distinct
        # genes -- and, with `require_core_role`, let a gene that is not
        # really there satisfy the core requirement.
        own_hits = [
            h for h in cluster.hits
            if h.family_key == family.key and h.superseded_by is None
        ]
        distinct_genes = {h.gene_name for h in own_hits}
        if len(distinct_genes) < floor.min_hits:
            continue
        if floor.require_core_role:
            # SearchHit.role already carries "core_MAT | flanking_conserved |
            # flanking_variable" directly (search.py) -- filter on the HIT's
            # own role, not on whether the gene NAME happens to be one the
            # family defines as core_MAT, which would not actually test
            # anything (a hit's gene_name is only ever one the family
            # declares in the first place).
            own_hits = [h for h in own_hits if h.role == "core_MAT"]
            if not own_hits:
                continue
        if floor.min_identity is not None and max(h.identity for h in own_hits) < floor.min_identity:
            continue
        admitted.append(family)
    return admitted


def _append_diagnostics_row(out_path: Path, row: dict) -> None:
    """Append one JSON line, never raising.

    Diagnostics are best-effort: losing the corpus is bad, losing a detection
    run because a log path was wrong is worse. The file is opened in APPEND
    mode and never truncated, which is what a batch writing many genomes to
    one file needs -- and why every row carries a `run_id`, so an accidental
    second run into the same file is detectable and de-duplicable instead of
    silently doubling the corpus.
    """
    try:
        with out_path.open("a") as f:
            f.write(json.dumps(row) + "\n")
    except OSError as exc:
        key = str(out_path)
        if key not in _DIAGNOSTICS_WRITE_FAILURES_LOGGED:
            _DIAGNOSTICS_WRITE_FAILURES_LOGGED.add(key)
            logger.warning("could not write evidence diagnostics to %s: %s", out_path, exc)


def _relaxed_results(
    clusters: list[GeneCluster],
    families: list[Family],
    searchable_genes: dict[FamilyKey, set[str]],
    evidence_floor: EvidenceFloor,
    ambiguity_floor: float = 0.5,
    polish_by: dict[tuple[int, FamilyKey, str], PolishOutcome] | None = None,
    contig_lengths: dict[str, int] | None = None,
) -> list[DetectionResult]:
    """Sub-floor clusters admitted on gene COUNT rather than gene fraction.

    The curator's relaxed second pass. A genuinely fragmented locus -- one
    split by a contig break, so its flanking genes sit on other contigs --
    can fall below the fraction floor through no fault of its own, while
    still showing a core gene and a real partner at high identity. The bar is
    "a core gene plus another gene, not necessarily a flank", which is
    exactly `EvidenceFloor(min_hits=2, require_core_role=True)`: one
    curator-ruled bar serving both polish admission and relaxed reporting, so
    there is a single number to defend rather than a second magic constant.

    Because it is a COUNT and not a fraction, a lone gene never passes however
    good it looks -- which is what keeps out the 22 sweep genera whose only
    evidence was a single HMG gene matched by both `sexM` and `sexP`.

    The caller runs this ONLY when the strict pass produced nothing
    genome-wide, so a genome with any confident call is untouched. Clusters
    that would clear the fraction floor are skipped anyway, since reporting
    one here would duplicate a strict result.

    Scope, measured: for Mucoromycota this admits nothing new. With
    searchable-only denominators (Plus 4, Minus 3) a core gene plus one other
    already scores 0.50 or 0.667 and clears the strict floor. It earns its
    place in phyla with richer rosters, where a real core+flank pair scores
    2/8 = 0.25 and is rejected today.
    """
    results: list[DetectionResult] = []
    families_by_key = {f.key: f for f in families}
    for cluster in clusters:
        admitted = _families_meeting_evidence_floor(cluster, families, evidence_floor)
        if not admitted:
            continue
        scores = score_cluster(cluster, families, searchable_genes=searchable_genes)
        admitted_keys = {f.key for f in admitted}
        for score in scores:
            if score.family_key not in admitted_keys:
                continue
            if score.fraction_found >= ambiguity_floor:
                continue  # the strict pass's own business, not this one's
            family = families_by_key[score.family_key]
            # Per-gene evidence is NOT optional. A locus reported with only a
            # list of gene names cannot be checked by anyone: in the 44-genus
            # sweep, 81 relaxed calls listed sexP and sexM with no coordinates
            # behind them, which made a genuine homothallic locus -- both
            # idiomorphs present, as documented for Syzygites -- impossible to
            # tell from two unrelated spurious HMG hits tens of kb apart.
            evidence = _gene_evidence([cluster], score.family_key, polish_by or {})
            segments = _segments_for([cluster], contig_lengths or {}, evidence)
            results.append(DetectionResult(
                family_key=score.family_key,
                contig=segments[0].contig,
                start=segments[0].start,
                end=segments[0].end,
                # Capped, never high: this call failed the fraction floor. Not
                # forced to low either -- a core gene plus a conserved flank at
                # high identity next to a contig break is real evidence, and
                # flattening it to low would conflate weak evidence with a
                # fragmented assembly.
                confidence="medium",
                idiomorph=assign_idiomorph(family, score.genes_found, _own_live_hits(cluster, family.key)),
                idiomorph_candidates=idiomorph_candidates(
                    family, score.genes_found, _own_live_hits(cluster, family.key)
                ),
                ambiguous_with=[],
                genes_found=score.genes_found,
                genes_missing=score.genes_missing,
                fragmented=False,
                genes_not_searchable=score.genes_not_searchable,
                detection_pass="relaxed",
                # A relaxed call never claims `mat_locus` -- curator's ruling,
                # 2026-09-21. `detection_pass` above still records the route,
                # so the two populations stay distinguishable.
                locus_class=apply_partial_locus(
                    classify_locus(cluster, family),
                    fraction_found=score.fraction_found,
                    ambiguity_floor=ambiguity_floor,
                    relaxed=True,
                ),
                span_exceeds_plausible_bound=span_exceeds_plausible_bound(
                    segments[0].start, segments[0].end, family
                ),
                segments=segments,
                gene_evidence=evidence,
                reference_records=sorted({e.reference_record_id for e in evidence}),
            ))
    return results


def _rescue_genes_in_idiomorph_scope(
    cluster: GeneCluster, family: Family, rescue_genes: set[str]
) -> tuple[set[str], set[str]]:
    """Split `rescue_genes` into (keep, skip) by the idiomorph this cluster
    already evidences. Returns `(rescue_genes, set())` whenever it cannot
    narrow safely.

    Rescue polish is the dominant cost in a genome-only run: for every
    (cluster, admitted family) the loop windows-polishes EVERY core_MAT gene
    still missing, with both exonerate and miniprot, at ~0.9 s a pair. On a
    two-idiomorph roster roughly half of those genes belong to the idiomorph
    this cluster is not. A MAT1-1 cluster cannot also hold MAT1-2-1, so
    polishing for it is guaranteed-futile work.

    Narrowing happens ONLY when the cluster's live, informative hits evidence
    exactly ONE idiomorph. Zero (nothing informative yet, or flanking-only
    evidence -- the localize-by-flanks case this project relies on) and two or
    more (a real homothallic both-idiomorphs locus, normal in curated records)
    both fall through unchanged. Genes carrying no `present_in_idiomorphs` are
    idiomorph-agnostic and are never skipped.

    NOT loss-free in principle, and the exposure is specific: the cluster's
    evidenced idiomorph could itself be wrong. The known way that happens is
    the shared-HMG cross-match (MAT1-1-3 vs MAT1-2-1, opposite idiomorphs,
    E 2.5e-17 to each other), where a single mis-attributed hit would then
    suppress the rescue that could have corrected it. Superseded hits are
    already excluded, which removes the cases the overlap resolver has
    caught; a cross-match that never overlapped is not covered. Skips are
    therefore recorded to the evidence-diagnostics stream so the real rate,
    and any call that changes because of them, can be measured rather than
    assumed.
    """
    evidenced = evidenced_idiomorphs(family, _own_live_hits(cluster, family.key))
    if len(evidenced) != 1:
        return rescue_genes, set()
    idiomorph = next(iter(evidenced))
    keep, skip = set(), set()
    for gene in family.genes:
        if gene["name"] not in rescue_genes:
            continue
        declared = gene.get("present_in_idiomorphs") or ()
        if declared and idiomorph not in declared:
            skip.add(gene["name"])
        else:
            keep.add(gene["name"])
    # A rescue gene the roster does not define cannot happen (`_missing_core_genes`
    # reads the roster), but keep the set total rather than silently shrinking.
    keep |= rescue_genes - keep - skip
    return keep, skip


def _write_polish_scope_diagnostics(
    out_path: Path,
    cluster: GeneCluster,
    family: Family,
    idiomorph: str,
    skipped: set[str],
    attempted: int,
    run_id: str,
    genome_id: str,
) -> None:
    """One JSON line per (cluster, family) whose rescue set was narrowed.

    `attempted` is the number of polish pairs actually run for this
    (cluster, family) after narrowing, so saved/attempted is computable per
    genome without re-deriving it from wall-clock time.
    """
    _append_diagnostics_row(out_path, {
        "kind": "polish_scope",
        "run_id": run_id,
        "genome_id": genome_id,
        "family": f"{family.key.phylum}:{family.key.locus_name}",
        "contig": cluster.contig, "cluster_start": cluster.start, "cluster_end": cluster.end,
        "evidenced_idiomorph": idiomorph,
        "rescues_skipped": sorted(skipped),
        "polish_pairs_attempted": attempted,
    })


def _write_evidence_diagnostics(
    out_path: Path,
    cluster: GeneCluster,
    family: Family,
    admitted: bool,
    run_id: str,
    genome_id: str,
) -> None:
    """Append one JSON line describing this (cluster, family) admission decision.

    The real calibration dataset `EvidenceFloor`'s docstring refers to. It
    could not actually serve that purpose before carrying `genome_id`: a row
    named the family, contig and cluster span but not the organism, so
    concatenating a batch's rows left no way to compute any per-genome
    statistic from them.
    """
    own_hits = [h for h in cluster.hits if h.family_key == family.key]
    _append_diagnostics_row(out_path, {
        "kind": "evidence",
        "run_id": run_id,
        "genome_id": genome_id,
        "family": f"{family.key.phylum}:{family.key.locus_name}",
        "contig": cluster.contig, "cluster_start": cluster.start, "cluster_end": cluster.end,
        "gene_count": len({h.gene_name for h in own_hits}),
        "hit_count": len(own_hits),
        "roles": sorted({h.role for h in own_hits}),
        "best_identity": max((h.identity for h in own_hits), default=None),
        "admitted": admitted,
    })


def _write_idiomorph_diagnostics(
    out_path: Path,
    resolution: IdiomorphResolution,
    family: Family,
    run_id: str,
    genome_id: str,
) -> None:
    """Append one JSON line per collapsed idiomorph pair.

    These rows are what the provisional overlap threshold
    (`idiomorph.DEFAULT_MIN_OVERLAP_FRACTION`, 0.5) and the provisional margin
    (`family_registry.DEFAULT_MIN_IDIOMORPH_MARGIN`, 5.0) are meant to be
    revised from, so both members' identity AND coverage are recorded, not
    just the verdict. Coverage in particular is the open question: the
    artifact being resolved is a shared protein domain, so a cross-hit should
    cover only part of its reference while the true gene covers all of it,
    which makes coverage the biologically motivated discriminator and identity
    a proxy that happened to score 23/23. It cannot be the rule today because
    it is `None` on the tblastn and exonerate paths -- which is itself
    something this corpus will show the rate of.
    """
    _append_diagnostics_row(out_path, {
        "kind": "idiomorph_resolution",
        "run_id": run_id,
        "genome_id": genome_id,
        "family": f"{family.key.phylum}:{family.key.locus_name}",
        "contig": resolution.contig,
        "winner": resolution.winner,
        "loser": resolution.loser,
        "winner_identity": resolution.winner_identity,
        "loser_identity": resolution.loser_identity,
        "margin": resolution.margin,
        "overlap_fraction": resolution.overlap_fraction,
        "winner_coverage": resolution.winner_coverage,
        "loser_coverage": resolution.loser_coverage,
    })


@dataclass(frozen=True)
class _RescueScope:
    """Which families the one batched genome-wide `tblastn` rescue runs for, and
    which of the hits it returns may be kept.

    The `tblastn` call itself always queries a family's FULL expected gene set,
    so this object only decides post-hoc which returned hits are folded back in.
    """

    #: The families passed to `search_localize` (one batched call for all).
    families: list[Family]
    #: family -> genes this rescue may contribute for that family. `None` means
    #: "no gene restriction": the family had no fast-path cluster at all, so its
    #: whole expected gene set is in scope.
    genes_by_family: dict[FamilyKey, set[str] | None]
    #: family -> [(fast-path cluster, that cluster's own found gene names for
    #: THIS family)], used to keep a rescue hit from being grafted onto a
    #: cluster that already has that gene.
    clusters_by_family: dict[FamilyKey, list[tuple[GeneCluster, set[str]]]]

    def accepts(self, hit: SearchHit, max_gap: int) -> bool:
        """True when a batched-rescue hit is one this rescue was actually run for.

        Three independent conditions, all defence in depth against the batched
        call returning more than any one family asked for:

        1. the hit's family must be one the rescue ran for at all;
        2. its gene must be one that family was missing SOMEWHERE (the whole
           gene set when the family had no cluster);
        3. it must not land on/next to a fast-path cluster of that family that
           ALREADY has that gene. Condition 2 is now per-cluster, so gene X can
           be in scope purely because cluster B lacks it; without condition 3 a
           genome-wide hit for X sitting on cluster A -- which already found X
           -- would be merged into cluster A, making A's already-annotated gene
           a polish candidate and letting a failed polish there cap A's tier.
           Proximity is judged with the same `max_gap` `cluster_hits` uses, so
           the test is "would this hit merge into that cluster".

        KNOWN LIMITATION (contained elsewhere, not fixed here): condition 3 is
        a DIRECT-overlap test against the cluster's span as it stands now,
        while `cluster_hits` merges by CHAINING. Two rescue hits that each pass
        this test independently can, together, chain a hit into a cluster it
        does not directly overlap. Making this test chain-aware would require
        running it AFTER the definitive clustering pass, which is a larger
        restructuring. The damaging consequences are instead contained in the
        polish loop, which refuses to treat a gene the cluster already has an
        annotated hit for as a polish candidate. The residual is that a
        chained-in spurious HSP still joins the cluster's hit list and so
        widens its span and any polish window computed from it.
        """
        if hit.family_key not in self.genes_by_family:
            return False
        allowed = self.genes_by_family[hit.family_key]
        if allowed is not None and hit.gene_name not in allowed:
            return False
        for cluster, found in self.clusters_by_family.get(hit.family_key, ()):
            if hit.gene_name not in found or hit.contig != cluster.contig:
                continue
            if hit.start <= cluster.end + max_gap and hit.end >= cluster.start - max_gap:
                return False
        return True


def _localization_rescue_targets(
    fast_path_clusters: list[GeneCluster], families: list[Family]
) -> _RescueScope:
    """Which families need the batched genome-wide `tblastn` rescue, and for
    which of their genes.

    Two shapes qualify, and they are served by the SAME batched
    `search_localize` call (one `makeblastdb` + one `tblastn` per run, never
    one per family):

    * **Zero fast-path hits** -- no foothold at all, so there is no window to
      polish against. The whole family's gene set is in scope, recorded as
      `None` ("no gene restriction").
    * **Partial foothold** -- the proteome placed some of the family's genes
      but not all of its `core_MAT` ones. Only those specific missing genes
      are in scope. Without this, a family with PARTIAL annotation coverage
      would get LESS search than a family with NONE: the only thing looked at
      for its missing gene would be a narrow `+-(6*aa + 2000) bp` window
      around its existing cluster (~+-3kb for a 100-aa reference protein).
      That is exactly this project's own motivating case -- a pheromone
      receptor correctly annotated while the short pheromone precursor beside
      it is absent from the same genome's own annotation -- and the missing
      gene is not guaranteed to sit inside that narrow window.

    "Missing" is asked per CLUSTER, not per family: a family's genes are
    collected per fast-path cluster (`_missing_core_genes`) and the family's
    rescue gene set is the UNION of what each of its own clusters lacks. A
    family legitimately has more than one independent locus in one genome
    (tetrapolar unlinked loci; homothallic/heterothallic switching cassettes),
    so gene X being present in cluster A says nothing about cluster B, and the
    family-wide question "is X missing anywhere?" would have answered "no" and
    denied cluster B the genome-wide look. Which cluster a kept hit ends up in
    is then decided spatially by `cluster_hits`, not by which cluster made the
    gene eligible.

    The narrow windowed rescue around the existing cluster still runs as
    before (see `_missing_core_genes` in the polish loop); this adds a
    genome-wide look for the same gene, whose hits flow into the ordinary
    clustering/polish machinery like any other localized hit.

    Gene scope is computed per family from that family's OWN fast-path hits,
    never pooled across families -- gene names are reused across families, so
    one family's found gene must never mask another's absence.
    """
    clusters_by_family: dict[FamilyKey, list[tuple[GeneCluster, set[str]]]] = {}
    for cluster in fast_path_clusters:
        for family in families:
            found = {h.gene_name for h in cluster.hits if h.family_key == family.key}
            if found:
                clusters_by_family.setdefault(family.key, []).append((cluster, found))

    rescue_families: list[Family] = []
    genes_by_family: dict[FamilyKey, set[str] | None] = {}
    for family in families:
        own = clusters_by_family.get(family.key)
        if not own:
            rescue_families.append(family)
            genes_by_family[family.key] = None  # no foothold: the whole family
            continue
        missing_core: set[str] = set()
        for cluster, _found in own:
            missing_core |= _missing_core_genes(cluster, family)
        if missing_core:
            rescue_families.append(family)
            genes_by_family[family.key] = missing_core
    return _RescueScope(rescue_families, genes_by_family, clusters_by_family)


def _contig_lengths(genome_fasta: Path) -> dict[str, int]:
    """Contig name -> length, or {} when the genome FASTA cannot be read.

    Used only to populate `contig_edge_distance` on fragmented multi-segment
    calls; an unreadable or absent FASTA degrades that field to None rather
    than failing the run.
    """
    try:
        from Bio import SeqIO

        return {record.id: len(record.seq) for record in SeqIO.parse(str(genome_fasta), "fasta")}
    except Exception:
        return {}


def span_exceeds_plausible_bound(start: int, end: int, family: Family) -> bool:
    """Is this locus wider than the family says a locus of its kind gets?

    A FLAG, not a filter, and public on purpose so a caller can ask the
    question without re-deriving the arithmetic. The curator ruled on
    2026-09-20 that wide loci are real and must still be found, so nothing in
    this pipeline may drop a call on this basis -- see
    `DEFAULT_MAX_PLAUSIBLE_LOCUS_SPAN_BP` for the measurements and for why the
    bound is 200 kb rather than the 120 kb first discussed.

    The comparison is strict (`>`): a locus exactly at the bound is within
    what the curator called acceptable.
    """
    return (end - start) > family.max_plausible_locus_span_bp


def _fragmented_family_segments(
    clusters: list[GeneCluster], family: Family
) -> list[GeneCluster] | None:
    """Clusters that together make one fragmented multi-contig locus, or None.

    A family is treated as fragmented across contigs only when all of:
      1. it has hits in clusters on more than one contig;
      2. NO single cluster contains all of its expected `core_MAT` genes
         (if one does, the other clusters are a second real locus copy --
         gene duplication and multi-allele co-occurrence are normal at MAT
         loci and must not be mislabelled as assembly fragmentation);
      3. the union of the chosen clusters DOES contain all of them, with each
         chosen cluster contributing at least one core gene the others lack.

    Condition 3 is a CONTRIBUTION test, not a disjointness test: the greedy
    cover below admits a cluster for contributing >=1 not-yet-covered core gene
    and never rejects it for ALSO holding a gene an already-chosen cluster has.
    Two chosen clusters can therefore each hold a real hit for the same gene
    name -- which is expected rather than exotic at MAT loci, where gene
    duplication and multi-allele co-occurrence are normal. That overlap is not
    an error here (each segment is a distinct genomic location and each hit is
    real), so `_gene_evidence` reports evidence per `(cluster, gene_name)`
    rather than collapsing the segments by gene name; see its docstring.
    """
    core_genes = {g["name"] for g in family.genes if g["role"] == "core_MAT"}
    if len(core_genes) < 2:
        return None

    own = [c for c in clusters if any(h.family_key == family.key for h in c.hits)]
    if len({c.contig for c in own}) < 2:
        return None

    found_in = {
        id(c): {h.gene_name for h in c.hits if h.family_key == family.key} & core_genes
        for c in own
    }
    if any(found == core_genes for found in found_in.values()):
        return None  # a complete single-contig call exists; this is not fragmentation

    chosen: list[GeneCluster] = []
    covered: set[str] = set()
    for cluster in sorted(own, key=lambda c: -len(found_in[id(c)])):
        new = found_in[id(cluster)] - covered
        if not new:
            continue
        chosen.append(cluster)
        covered |= new

    if covered != core_genes or len({c.contig for c in chosen}) < 2:
        return None
    return chosen


#: Default multiple of a curated reference protein's own nucleotide-equivalent
#: length (aa * 3) used as polish-window padding on each side of a cluster.
#: 2.0 is chosen so a gene up to three times its reference's coding length still
#: fits inside the window: `tblastn` frequently anchors on only a short conserved
#: domain of a divergent MAT protein (an HMG box or homeodomain), so a window
#: sized to the HSP alone would truncate the real gene.
DEFAULT_WINDOW_PROTEIN_LENGTH_MULTIPLE = 2.0
#: Default additive allowance for intronic sequence the reference protein's
#: length cannot account for. 2000 bp is generous for fungi, where the great
#: majority of introns are well under 500 bp; it also serves as the whole padding
#: when no curated length is known for a gene.
DEFAULT_WINDOW_MAX_INTRON_BP = 2_000


def _window_padding(
    protein_aa: int | None,
    protein_length_multiple: float = DEFAULT_WINDOW_PROTEIN_LENGTH_MULTIPLE,
    max_intron_bp: int = DEFAULT_WINDOW_MAX_INTRON_BP,
) -> int:
    """Padding, in bp, to add on each side of a cluster before polishing one gene.

    Derived from that gene's OWN curated reference protein length rather than a
    single hardcoded constant, per the spec: `protein_aa * 3 *
    protein_length_multiple + max_intron_bp`. Both knobs are named and
    overridable (see `run_pipeline`'s `window_*` parameters). `protein_aa` is
    None when the curated database carries no protein for that
    `(family, gene)`; the padding then degrades to `max_intron_bp` alone rather
    than to zero, so an unknown length never produces a window too tight to
    align in.
    """
    return int((protein_aa or 0) * 3 * protein_length_multiple) + max_intron_bp


def _padded_window(
    cluster: GeneCluster, padding: int, contig_lengths: dict[str, int]
) -> tuple[str, int, int]:
    """(contig, start, end) 1-based inclusive, clamped to the contig when its
    length is known (an unreadable genome FASTA leaves `contig_lengths` empty,
    in which case the un-clamped end is passed through -- `_extract_window`'s
    slice tolerates an end past the contig)."""
    start = max(1, cluster.start - padding)
    end = cluster.end + padding
    length = contig_lengths.get(cluster.contig)
    if length:
        end = min(end, length)
    return (cluster.contig, start, end)


def _own_model(
    model: PolishModel | None, family_key: FamilyKey, gene_name: str, contig: str
) -> PolishModel | None:
    """A polished model, or None when it is not this family's own named gene on
    the contig the polish window was opened on.

    Both polishing wrappers already drop a model whose gene name or curated
    record does not match what was asked for, so this is defence in depth --
    kept because the retired windowed second pass carried the same filter, and
    because a model silently attributed to the wrong family is precisely the
    cross-family bug class this file has shipped before. The `contig` check
    closes the same hole in the spatial dimension: a model whose coordinates
    name a different contig than the window it was asked for cannot belong to
    this cluster, and accepting it would write coordinates from somewhere else
    in the genome into this cluster's evidence. Dropping a model degrades the
    gene to `unpolished` (its raw localization hit still stands) rather than
    writing another family's or another contig's coordinates into this family's
    evidence.
    """
    if model is None:
        return None
    if model.family_key != family_key or model.gene_name != gene_name:
        return None
    if model.contig != contig:
        return None
    return model


def _hit_from_model(model: PolishModel) -> SearchHit:
    """A polished gene model re-expressed as a SearchHit, so a fast-path rescue's
    newly found gene can be fed back into its own cluster's hit list and counted
    by the existing `score_cluster` / fragmentation logic unchanged."""
    return SearchHit(
        family_key=model.family_key, gene_name=model.gene_name, role=model.role,
        contig=model.contig, start=model.start, end=model.end, strand=model.strand,
        identity=model.identity, reference_record_id=model.reference_record_id,
        method=model.method, coverage=None,
    )


def _model_dict(model: PolishModel) -> dict:
    """A `PolishModel` re-expressed as a plain, YAML/GFF3-attribute-friendly
    dict -- used only to carry the non-canonical tool's model on a
    `polished_disagree` `GeneEvidence` so a human reviewer can inspect the
    disagreement, never for the canonical fields (those stay as
    `GeneEvidence`'s own top-level attributes)."""
    return {
        "contig": model.contig,
        "start": model.start,
        "end": model.end,
        "strand": model.strand,
        "exons": [{"start": e.start, "end": e.end} for e in model.exons],
        "identity": model.identity,
        "method": model.method,
    }


#: Preference order for RAW (unpolished) hits of the same gene, most preferred
#: first. `diamond_proteome` is a protein-vs-protein match against a real
#: annotated gene model; `tblastn_genome` is a deliberately approximate,
#: localization-only alignment with no splice awareness. This is a named
#: method preference, NOT a score comparison: their identity numbers come from
#: different searches and are not on a common scale.
_RAW_METHOD_PREFERENCE = ("diamond_proteome", "tblastn_genome")


def _raw_method_rank(method: str) -> int:
    try:
        return _RAW_METHOD_PREFERENCE.index(method)
    except ValueError:
        return len(_RAW_METHOD_PREFERENCE)


def _live_hits_for_evidence(hits: list[SearchHit]) -> list[SearchHit]:
    """Drop superseded hits, unless that would leave a gene with no evidence.

    A gene counts as found because of a LIVE hit, so the evidence reported for
    it must be that hit. Measured on Syzygites sp. MES_3091 scaffold_11: sexP
    counted as found because of a hit at 148,515, but the report displayed the
    superseded overlapping hit at 143,389 instead -- the very hit that had been
    ruled out. The displayed evidence contradicted the scoring that admitted it.

    Per gene, not globally: a gene whose hits are ALL superseded keeps them, so
    the ambiguity stays visible in the report rather than vanishing. That gene
    is not in `genes_found` anyway, so this only affects what a reader sees.
    """
    by_gene: dict[str, list[SearchHit]] = {}
    for hit in hits:
        by_gene.setdefault(hit.gene_name, []).append(hit)
    kept: list[SearchHit] = []
    for gene_hits in by_gene.values():
        live = [h for h in gene_hits if h.superseded_by is None]
        kept.extend(live or gene_hits)
    return kept


def _gene_evidence(
    member_clusters: list[GeneCluster],
    family_key: FamilyKey,
    polish_by: dict[tuple[int, FamilyKey, str], PolishOutcome],
) -> list[GeneEvidence]:
    """Report-ready per-gene evidence for ONE family across its own clusters.

    Per gene, the canonical `PolishOutcome.canonical` model wins when that gene
    was polished in that cluster; a `STATUS_UNPOLISHED` gene (canonical is None)
    and any gene that was never a polish candidate at all -- e.g. a gene the
    diamond fast path already found -- fall back to that gene's best raw
    `SearchHit` (chosen by `_RAW_METHOD_PREFERENCE`, with identity breaking
    ties only within one method). That fallback is the only path by which
    a `GeneEvidence` is built from a `SearchHit` rather than a `PolishModel`.
    The two situations are reported with DIFFERENT statuses so a curator can
    tell them apart: `STATUS_UNPOLISHED` when a `PolishOutcome` exists for the
    gene (it WAS localized and sent through both polish tools, but neither
    confirmed it -- genuinely uncertain, caps the family's tier), versus
    `STATUS_NOT_POLISH_CANDIDATE` when no `PolishOutcome` was ever computed for
    it (it never entered the localize/polish pipeline at all -- solid, no
    uncertainty implied). This distinction is purely for the human-readable
    report: `_any_gene_unpolished` (tiering) reads `PolishOutcome.status`
    directly, never `GeneEvidence.status`, and is completely unaffected by it.

    `polish_by` is keyed `(id(cluster), family_key, gene_name)` and is read here
    ONLY for this family's own genes in these exact clusters, so a polish result
    from another family, or from an unrelated cluster of this same family, can
    never be attributed to this call.

    Selection is likewise per `(cluster, gene_name)`, NEVER per bare
    `gene_name`. `member_clusters` holds more than one cluster only for a
    fragmented multi-segment call, and those clusters are by construction
    distinct genomic locations. `_fragmented_family_segments`' greedy cover is a
    CONTRIBUTION test, not a disjointness test (a cluster is admitted for
    contributing >=1 not-yet-covered core gene; nothing rejects it for also
    sharing a gene name with an already-chosen cluster), and gene duplication /
    multi-allele co-occurrence is normal at MAT loci -- so two segments of one
    fragmented call can each hold a real, correct hit for the same gene name.
    Collapsing those into one entry by gene name silently discarded one
    segment's real gene model from the curation-facing report and left that
    segment's reported span wider than any gene listed for it. Each cluster
    therefore contributes its own evidence, built from its OWN hits and its OWN
    `PolishOutcome`s -- exactly as the single-cluster case has always done --
    and the per-segment dictionaries are concatenated, with no cross-cluster
    identity or polished-vs-raw comparison at any point. This covers every gene
    the fragmentation cover assigned to a segment (it is that segment's own hit
    that is reported) and additionally keeps a gene the cover did not assign to
    that segment, including non-core genes, which the cover never considers at
    all. Within ONE cluster the ordering is unchanged: a polished
    `PolishOutcome.canonical` wins outright over that cluster's raw hits, and
    two raw hits are ranked by `_RAW_METHOD_PREFERENCE` with identity breaking
    ties only inside one method, so identity is still never compared across
    tools.
    """
    best: dict[tuple[int, str], GeneEvidence] = {}
    for cluster in member_clusters:
        raw_by_gene: dict[str, SearchHit] = {}
        for hit in _live_hits_for_evidence(cluster.hits):
            if hit.family_key != family_key:
                continue
            current = raw_by_gene.get(hit.gene_name)
            # Ranking among THIS cluster's own raw SearchHits for one gene: the
            # named method preference decides first, and identity only breaks a tie
            # between two hits from the SAME method. A diamond `pident` and a
            # tblastn `pident` describe different alignments and must not be
            # ranked against each other as if they were one scale.
            if current is None:
                raw_by_gene[hit.gene_name] = hit
                continue
            hit_rank = _raw_method_rank(hit.method)
            current_rank = _raw_method_rank(current.method)
            if hit_rank != current_rank:
                if hit_rank < current_rank:
                    raw_by_gene[hit.gene_name] = hit
            elif hit.method == current.method and hit.identity > current.identity:
                raw_by_gene[hit.gene_name] = hit

        outcomes = {
            gene_name: outcome
            for (cluster_id, key, gene_name), outcome in polish_by.items()
            if cluster_id == id(cluster) and key == family_key
        }

        for gene_name in set(raw_by_gene) | set(outcomes):
            outcome = outcomes.get(gene_name)
            if outcome is not None and outcome.canonical is not None:
                model = outcome.canonical
                # The non-canonical tool's model is surfaced only when the two
                # tools actually disagreed -- a `polished_agree` pair's second
                # model is redundant (within tolerance of the canonical one)
                # and `polished_single` has no second model to show.
                alternate = None
                if outcome.status == STATUS_DISAGREE:
                    other = (
                        outcome.miniprot_model
                        if outcome.canonical is outcome.exonerate_model
                        else outcome.exonerate_model
                    )
                    alternate = _model_dict(other) if other is not None else None
                evidence = GeneEvidence(
                    gene_name=model.gene_name, role=model.role, contig=model.contig,
                    start=model.start, end=model.end, strand=model.strand,
                    identity=model.identity, coverage=None,
                    reference_record_id=model.reference_record_id, method=model.method,
                    status=outcome.status, alternate_model=alternate,
                    exons=tuple((e.start, e.end) for e in model.exons) if model.exons else None,
                )
            else:
                hit = raw_by_gene.get(gene_name)
                if hit is None:
                    continue  # unpolished with nothing localized: no evidence to report
                # `outcome` is None exactly when this gene never entered the
                # localize/polish pipeline at all (e.g. found only via the
                # fast-path diamond hit, or an already-present core gene) --
                # a solid result with no implied uncertainty, tagged
                # STATUS_NOT_POLISH_CANDIDATE. `outcome is not None` here
                # means it WAS a polish candidate (localized/rescued and sent
                # through both tools) but neither tool produced a model
                # (outcome.canonical is None) -- genuinely uncertain, tagged
                # STATUS_UNPOLISHED, exactly what caps the family's tier via
                # `_any_gene_unpolished` (which reads PolishOutcome.status
                # directly and is unaffected by this GeneEvidence.status
                # split). In both cases the raw SearchHit stands as this
                # gene's evidence with no polished model behind it.
                status = STATUS_UNPOLISHED if outcome is not None else STATUS_NOT_POLISH_CANDIDATE
                evidence = GeneEvidence(
                    gene_name=hit.gene_name, role=hit.role, contig=hit.contig,
                    start=hit.start, end=hit.end, strand=hit.strand,
                    identity=hit.identity, coverage=hit.coverage,
                    reference_record_id=hit.reference_record_id, method=hit.method,
                    status=status, alternate_model=None,
                )
            # Keyed per (cluster, gene): a cluster visits each of its own gene
            # names exactly once, so this never overwrites and no evidence from
            # another segment can displace this one.
            best[(id(cluster), gene_name)] = evidence
    # `gene_name` is part of the sort key so two entries that happen to share a
    # contig and a start coordinate still order deterministically.
    return sorted(best.values(), key=lambda e: (e.contig, e.start, e.end, e.gene_name))


def _segments_for(
    clusters: list[GeneCluster],
    contig_lengths: dict[str, int],
    evidence: list[GeneEvidence],
) -> list[LocusSegment]:
    """One `LocusSegment` per member cluster, each spanning the UNION of that
    cluster's own frozen span and every gene this same result reports on that
    cluster's contig.

    A `GeneCluster`'s span is frozen at `cluster_hits` time from the raw
    localization HSPs, but `_gene_evidence` reports polished/rescued coordinates
    that can extend past it -- a rescued gene at 150-260 inside a cluster whose
    frozen span is 300-400. Reporting the frozen span alone produces an invalid
    GFF3: a child `gene` feature outside its parent `MAT_locus` feature's
    declared range and outside the `##sequence-region` pragma. The spec requires
    the reported window to cover the union of all its genes' canonical
    coordinates rather than be silently truncated.

    Widening is done per contig, not per cluster: when one result reports two
    clusters on the SAME contig (possible for a fragmented call whose chosen
    clusters span >=2 contigs but include two on one of them) both are widened
    by the same contig's genes and can overlap. That is deliberate -- an
    over-wide reported span still contains every gene it claims, whereas a
    per-cluster nearest-gene assignment would have to guess which cluster a
    polished coordinate belongs to.
    """
    segments = []
    for cluster in sorted(clusters, key=lambda c: (c.contig, c.start)):
        on_contig = [e for e in evidence if e.contig == cluster.contig]
        start = min([cluster.start] + [e.start for e in on_contig])
        end = max([cluster.end] + [e.end for e in on_contig])
        length = contig_lengths.get(cluster.contig)
        edge = min(start - 1, length - end) if length else None
        segments.append(LocusSegment(cluster.contig, start, end, edge))
    return segments


def run_pipeline(
    genome_fasta: Path,
    proteome_fasta: Path | None,
    taxid: int | None,
    db_root: Path,
    reference_fasta: Path,
    search_fast_path: Callable = search_fast_path,
    search_localize: Callable = search_localize,
    polish_with_exonerate: Callable = polish_with_exonerate,
    polish_with_miniprot: Callable = polish_with_miniprot,
    max_gap: int | None = None,
    ambiguity_floor: float = 0.5,
    short_orf_aa_floor: int = 60,
    idiomorph_overlap_fraction: float = DEFAULT_MIN_OVERLAP_FRACTION,
    relaxed_second_pass: bool = True,
    window_protein_length_multiple: float = DEFAULT_WINDOW_PROTEIN_LENGTH_MULTIPLE,
    window_max_intron_bp: int = DEFAULT_WINDOW_MAX_INTRON_BP,
    polish_tolerance_bp: int = 10,
    evidence_floor: EvidenceFloor = EvidenceFloor(),
    evidence_diagnostics_path: Path | None = None,
    routing: RoutingDecision | None = None,
    allow_cross_contig_fragments: bool = False,
    #: How many of a locus's genes must rest on a real gene model for it to be
    #: reported. See `MIN_POLISHED_GENES` for the measurement behind the
    #: default and `_modelled_gene_count` for what counts. 0 disables the bar,
    #: which is what the tests written before it use to keep asserting the
    #: behaviour they were written for.
    min_polished_genes: int = MIN_POLISHED_GENES,
    #: Curated records withheld from THIS run, for leave-one-out recall. The
    #: caller must build `reference_fasta` with the same set (see
    #: `reference_fasta.build_reference_fasta`); this parameter additionally
    #: keeps them out of every db-derived quantity -- polish-window padding and
    #: the short-ORF/unsearchable split -- so a held-out answer cannot shape a
    #: run that is meant to be blind to it. `detect.holdout` chooses the set,
    #: and its RADIUS is what makes the resulting number meaningful.
    exclude_record_ids: frozenset[str] = frozenset(),
    #: NCBI translation table for THIS genome. None means "derive from the
    #: taxid"; an explicit value wins, exactly as `--phylum` overrides taxid
    #: routing. Falls back to 1 when the taxonomy lookup cannot answer -- never
    #: guessed. The CUG-Ser1 clade (Serinales) is table 12 and reads CTG as
    #: serine; translating those genomes with table 1 is systematically wrong.
    genetic_code: int | None = None,
    genetic_code_resolver: Callable[[int], int | None] = default_genetic_code,
) -> DetectionOutcome:
    # `routing` lets the caller route ONCE and reuse the decision, which the
    # CLI must do: it has to know the routed families BEFORE this call, so it
    # can restrict `build_reference_fasta`'s query set to them. Routing here as
    # well would repeat the taxonomy lookup and, worse, could disagree with the
    # reference FASTA the caller already built. When it is omitted the old
    # self-routing behaviour is unchanged.
    if routing is None:
        routing = route(taxid, load_all_families(db_root))
    families = routing.families
    # The clustering gap is per-locus CURATION data (`order.yml`
    # `max_cluster_gap_bp`), not a code constant: how spread out a MAT locus is
    # differs by clade, and the curator ruled on 2026-09-20 that 25 kb fits
    # Ascomycota but is too tight for Mucoromycota. The run's gap is the
    # MAXIMUM over the routed families because the two errors are not
    # symmetric: too LARGE a gap under-splits, which is recoverable -- the
    # evidence floor and polishing still discriminate gene by gene inside an
    # over-large cluster -- whereas too SMALL a gap over-splits, silently
    # destroying a real locus by cutting it in two, and nothing downstream can
    # rejoin the halves. An explicit `max_gap` from the caller still wins; None
    # means derive.
    #
    # The routing mode qualifies the maximum: on `phylum_fallback` and
    # `exhaustive` the default stands instead, because a failed route gives no
    # basis for inheriting the widest locus's gap. See
    # `derive_max_cluster_gap`'s docstring for the measured reason (the 120 kb
    # Tremellales MAT locus would otherwise cluster every unrouted
    # Basidiomycota genome at 120 kb).
    if max_gap is None:
        max_gap = derive_max_cluster_gap(families, routing_mode=routing.routing_mode)
    # Derived from the taxid when not supplied, out of the SAME cached efetch
    # document the router just read, so this costs no extra network call. A
    # failed lookup degrades to the standard table rather than raising.
    genetic_code_error = None
    if genetic_code is None and taxid is not None:
        try:
            genetic_code = genetic_code_resolver(taxid)
        except Exception as exc:
            genetic_code = None
            # Recorded, not just absorbed: a CTG-clade yeast (table 12)
            # translated with table 1 is a different search.
            genetic_code_error = f"{type(exc).__name__}: {str(exc)[:160]}; used table 1"
    if genetic_code is None:
        genetic_code = 1
    record_families = load_record_families(db_root)
    #: (family, gene) -> the roster's `polish` setting; absent means both tools.
    polish_mode = {(f.key, g["name"]): g.get("polish") for f in families for g in f.genes}
    record_protein_lengths: dict[tuple[str, str], int] = {}
    protein_lengths = _curated_protein_lengths(
        db_root, families, record_families, exclude_record_ids,
        by_record=record_protein_lengths)
    short_orf_by_family = _short_orf_genes(
        db_root, families, record_families, short_orf_aa_floor, exclude_record_ids
    )
    families_by_key = {f.key: f for f in families}

    # Stage 0/1 -- at most two searches per run, each batched over many families.
    # The genome-only path localizes with a single batched genome-wide tblastn
    # call covering every routed family. The fast path uses the supplied
    # proteome, then falls back to ONE batched tblastn localization covering
    # both every family the proteome produced no hit for at all (no cluster, so
    # no window to polish against) and every partially-annotated family's own
    # specific missing core_MAT genes. Without the latter, a family with PARTIAL
    # annotation coverage would get LESS search than one with none -- see
    # `_localization_rescue_targets`. Either way this is the blind spot (a gene
    # absent from a genome's own annotation) this pipeline exists to close.
    hits: list[SearchHit] = []
    # id() of every hit that came from tblastn localization (either path). Used
    # to decide, per cluster and family, which genes need their approximate HSP
    # coordinates refined. Keying on identity rather than on the hit's `method`
    # string means an injected/stubbed search cannot accidentally be classified
    # by what it happened to name its method. `hits` (and, after clustering,
    # `cluster.hits`) holds every one of these objects alive for the whole
    # function, so no id can be recycled.
    localized_hit_ids: set[int] = set()
    if proteome_fasta is not None:
        hits.extend(search_fast_path(proteome_fasta, families, reference_fasta, record_families))
        # Cluster the fast-path hits FIRST, purely to answer "which genes is
        # each individual existing locus of this family missing?". This is a
        # pure in-memory grouping (`cluster_hits`), not another search, and it
        # is what makes rescue eligibility per (family, gene, cluster) instead
        # of per (family, gene): a family with an unlinked second locus must
        # get a genome-wide look for a gene its OTHER locus already has. The
        # definitive clustering still happens once below, over the fast-path
        # and rescued hits together, so a rescued hit that lands inside an
        # existing cluster joins it and one that lands elsewhere seeds a new
        # cluster -- exactly as before this change.
        rescue_scope = _localization_rescue_targets(
            cluster_hits(hits, max_gap=max_gap), families
        )
        if rescue_scope.families:
            rescued = [
                h
                for h in search_localize(
                    genome_fasta, rescue_scope.families, reference_fasta,
                    record_families, genetic_code=genetic_code,
                )
                # Defence in depth: only the families this rescue was actually
                # run for may gain hits from it, only for the genes some cluster
                # of theirs was missing, and never on top of a cluster that
                # already has that gene (see `_RescueScope.accepts`).
                if rescue_scope.accepts(h, max_gap)
            ]
            localized_hit_ids.update(id(h) for h in rescued)
            hits.extend(rescued)
    else:
        localized = search_localize(
            genome_fasta, families, reference_fasta, record_families,
            genetic_code=genetic_code,
        )
        localized_hit_ids.update(id(h) for h in localized)
        hits.extend(localized)

    # Drop alignments too short to carry information BEFORE clustering, so a
    # fragment can never contribute a gene name to a cluster or inflate the
    # evidence floor's gene count. Measured need: a 27 bp "sexM" (nine codons)
    # and a 48/51 bp sexM/sexP pair were each being reported as loci.
    hits = drop_low_quality_hits(hits)

    clusters = cluster_hits(hits, max_gap=max_gap)

    # Resolve overlapping mutually-exclusive idiomorph genes HERE, before
    # anything counts distinct genes. sexM and sexP share an HMG box, so one
    # real locus gene draws both references; left unresolved it inflates the
    # evidence floor's gene count, stops `expected_genes_for_idiomorph`
    # narrowing the roster, and leaves the idiomorph uncallable. Deferring
    # this to scoring would fix the last two and leave the first.
    #
    # Per cluster and per family: two hits only describe the same gene if
    # they are in the same cluster, and gene names are only mutually
    # exclusive within the family that declares their idiomorphs.
    idiomorph_events_by_cluster: dict[int, list[IdiomorphResolution]] = {}
    #: (family, event) for every resolution in the run, for the diagnostics
    #: corpus. Paired at resolution time rather than re-derived afterwards:
    #: a cluster routinely mixes families, and gene names are not unique
    #: across them, so recovering the family from the event alone is guesswork.
    resolutions_with_family: list[tuple[Family, IdiomorphResolution]] = []
    resolved_clusters: list[GeneCluster] = []
    for cluster in clusters:
        cluster_hits_out = list(cluster.hits)
        events: list[IdiomorphResolution] = []
        for family in families:
            own = [h for h in cluster_hits_out if h.family_key == family.key]
            if len(own) < 2:
                continue
            resolved_own, own_events = resolve_idiomorph_overlaps(
                own, family, min_overlap_fraction=idiomorph_overlap_fraction
            )
            if not own_events:
                continue
            by_id = dict(zip((id(h) for h in own), resolved_own))
            cluster_hits_out = [by_id.get(id(h), h) for h in cluster_hits_out]
            events.extend(own_events)
            resolutions_with_family.extend((family, e) for e in own_events)
        new_cluster = GeneCluster(
            cluster.contig, cluster.start, cluster.end, cluster_hits_out
        )
        if events:
            idiomorph_events_by_cluster[id(new_cluster)] = events
        resolved_clusters.append(new_cluster)
    clusters = resolved_clusters

    # One id per run, so a second run appending into the same diagnostics file
    # is detectable rather than silently doubling the corpus, and one per
    # genome, so a batch's rows can be told apart at all.
    run_id = uuid.uuid4().hex[:12]
    genome_id = genome_fasta.stem
    if evidence_diagnostics_path is not None:
        for family, event in resolutions_with_family:
            _write_idiomorph_diagnostics(
                evidence_diagnostics_path, event, family,
                run_id=run_id, genome_id=genome_id,
            )

    # What this run could have found, read back from the FASTA it searched
    # with. Composed with the short-ORF exclusion because a gene that cannot
    # be found is a gene that cannot be found, whichever reason applies: both
    # must leave the denominator, or `fraction_found`'s ceiling tracks gaps in
    # the curated database instead of the biology.
    searchable_genes = {
        key: genes - short_orf_by_family.get(key, set())
        for key, genes in searchable_genes_by_family(
            reference_fasta, record_families
        ).items()
    }

    # Read contig lengths unconditionally so contig_edge_distance is populated
    # (or left None on an unreadable FASTA) consistently for every segment of
    # every result, independent of whether some OTHER family in this same run
    # happened to be fragmented. Gating this on `fragmented_segments` used to
    # make two otherwise-identical single-contig runs disagree on
    # contig_edge_distance purely because of an unrelated family elsewhere in
    # the genome. It is read before polishing because the polish window is
    # clamped to the contig's real end.
    contig_lengths = _contig_lengths(genome_fasta)

    # Stage 2/3 -- polish. `polish_by` is keyed by
    # (id(cluster), family.key, gene_name), NOT by family.key or gene name
    # alone. The same family can have several independent spatial clusters in
    # one genome (gene duplication / multi-allele co-occurrence is normal at MAT
    # loci) and the same gene name is reused across families, so a result must
    # be attributable to exactly one cluster, one family and one gene. This is
    # the same scoping the retired second-pass tracking set used, now carrying a
    # per-gene outcome instead of a per-family boolean.
    polish_by: dict[tuple[int, FamilyKey, str], PolishOutcome] = {}
    for cluster in clusters:
        admitted_families = _families_meeting_evidence_floor(cluster, families, evidence_floor)
        if evidence_diagnostics_path is not None:
            # Diagnostics cover every family the OLD unconditional-admit set
            # would have considered (i.e. every family with >=1 own hit in
            # this cluster), not just the ones the real `evidence_floor`
            # actually admits -- so both admitted and would-have-been-rejected
            # cases are captured for later threshold analysis.
            admitted_keys = {f.key for f in admitted_families}
            for family in _families_meeting_evidence_floor(
                cluster, families, _DIAGNOSTICS_CANDIDATE_FLOOR
            ):
                _write_evidence_diagnostics(
                    evidence_diagnostics_path, cluster, family,
                    admitted=family.key in admitted_keys,
                    run_id=run_id, genome_id=genome_id,
                )
        for family in admitted_families:
            # Eligibility is decided per (cluster, family), NOT from a single
            # global "is this a genome-only run" flag. Since the fast path can
            # now also carry tblastn-localized clusters (the zero-hit rescue
            # above), a global flag would leave a rescued cluster's genes
            # unpolished and would stop an unpolished gene there from ever
            # capping its family's tier.
            #
            # Genes this family ALREADY has in THIS cluster from a real,
            # non-localized hit (a fast-path diamond match against an annotated
            # gene model). Such a gene needs no refining here: its coordinates
            # come from an annotated model, not from an approximate HSP.
            #
            # Excluding them is the containment for the guard/chaining mismatch
            # in `_RescueScope.accepts`. `accepts` tests a rescued hit against a
            # cluster's ORIGINAL span (+-max_gap), but `cluster_hits` merges by
            # CHAINING: hit1-hit2 and hit2-hit3 each within max_gap puts hit1
            # and hit3 in one cluster even when they are much further apart. So
            # a spurious rescue hit for gene X, correctly judged by `accepts` as
            # NOT sitting on cluster A, can still be chained into cluster A by a
            # second, separately-accepted rescue hit that bridges the gap -- and
            # X was only in the family's rescue scope because some OTHER cluster
            # of the family lacked it. Treating that chained-in HSP as a
            # localized gene of A made it a polish candidate there; a failed
            # polish then capped A's tier via `_any_gene_unpolished`, and a
            # successful one would have overridden A's real annotated
            # coordinates in `_gene_evidence`, which ranks a cluster's polished
            # model above that same cluster's raw hits for the gene (a
            # WITHIN-cluster rule; across the segments of a fragmented call
            # evidence is no longer collapsed at all). Dropping such genes makes both effects
            # impossible for a gene the cluster already genuinely has.
            #
            # This never suppresses a real rescue: a gene genuinely missing from
            # this cluster has no non-localized hit here, so it is not in this
            # set, and the zero-hit and partial-foothold rescue paths are
            # untouched (on a genome-only run nothing is non-localized, so the
            # set is empty). `rescue_genes` below is likewise unaffected, since
            # a gene present in the cluster is by definition not missing from it.
            already_annotated_genes = {
                h.gene_name
                for h in cluster.hits
                if h.family_key == family.key and id(h) not in localized_hit_ids
            }
            # Localized genes: this family's own genes that tblastn placed in
            # THIS cluster, whose approximate HSP coordinates need refining.
            localized_genes = {
                h.gene_name
                for h in cluster.hits
                if h.family_key == family.key and id(h) in localized_hit_ids
            } - already_annotated_genes
            # Rescue genes: this family's own core_MAT genes still missing from
            # this cluster, checked strictly against this family's own hits in
            # this cluster (see _missing_core_genes), so one family's presence
            # never masks another family's absence. Nothing localized these, so
            # a failed rescue leaves the gene genuinely missing rather than
            # "unpolished" (see below).
            rescue_genes = _missing_core_genes(cluster, family) - localized_genes
            # Cut 1 of the polish-cost work: a rescue for a gene belonging to
            # the idiomorph this cluster is NOT cannot succeed, so do not pay
            # 0.9 s to find that out. Narrows only on unambiguous evidence;
            # see `_rescue_genes_in_idiomorph_scope` for the exposure.
            rescue_genes, skipped_rescues = _rescue_genes_in_idiomorph_scope(
                cluster, family, rescue_genes
            )

            if skipped_rescues and evidence_diagnostics_path is not None:
                _write_polish_scope_diagnostics(
                    evidence_diagnostics_path, cluster, family,
                    idiomorph=next(iter(
                        evidenced_idiomorphs(family, _own_live_hits(cluster, family.key))
                    )),
                    skipped=skipped_rescues,
                    attempted=len(localized_genes | rescue_genes),
                    run_id=run_id, genome_id=genome_id,
                )

            for gene_name in sorted(localized_genes | rescue_genes):
                rescue = gene_name not in localized_genes
                padding = _window_padding(
                    protein_lengths.get((family.key, gene_name)),
                    protein_length_multiple=window_protein_length_multiple,
                    max_intron_bp=window_max_intron_bp,
                )
                window = _padded_window(cluster, padding, contig_lengths)
                # Both tools get the SAME window and the same gene, so their
                # models are directly comparable.
                # A roster gene marked `polish: miniprot` skips exonerate:
                # curator's ruling, 2026-09-26, for the long PAP1/OBP1/PIK1
                # flanks, where exonerate was ~90% of runtime and a 60-genome
                # ablation gave identical genotypes on miniprot alone.
                # `classify` then records a `polished_single` model.
                if polish_mode.get((family.key, gene_name)) == "miniprot":
                    exonerate_model = None
                else:
                    exonerate_model = polish_with_exonerate(
                        genome_fasta=genome_fasta, family=family, gene_name=gene_name,
                        reference_fasta=reference_fasta, record_families=record_families,
                        window=window, genetic_code=genetic_code,
                    )
                miniprot_model = polish_with_miniprot(
                    genome_fasta=genome_fasta, family=family, gene_name=gene_name,
                    reference_fasta=reference_fasta, record_families=record_families,
                    window=window, genetic_code=genetic_code,
                )
                outcome = classify(
                    _own_model(exonerate_model, family.key, gene_name, window[0]),
                    _own_model(miniprot_model, family.key, gene_name, window[0]),
                    polish_tolerance_bp,
                )
                if rescue:
                    if outcome.canonical is None:
                        # A rescue that found nothing is not an "unpolished"
                        # gene: nothing localized it, so there is no rough hit
                        # to fall back on and the gene stays genuinely missing.
                        # Recording it as unpolished would wrongly cap the tier
                        # of a family whose gene was never evidenced at all.
                        continue
                    cluster.hits.append(_hit_from_model(outcome.canonical))
                polish_by[(id(cluster), family.key, gene_name)] = outcome

    # Resolve overlapping idiomorph pairs AGAIN, now that polishing has run.
    # The first pass (before polishing) is what keeps the evidence floor
    # honest, but polishing MOVES coordinates: exonerate/miniprot replace an
    # approximate tblastn span with a refined gene model, and a pair that did
    # not overlap beforehand can overlap afterwards. Measured on Syzygites sp.
    # MES_3091, three of five reported loci had post-polish sexM/sexP overlaps
    # of 100%, 32 bp and 47 bp that the pre-polish pass never saw, so each was
    # reported as two genes and left idiomorph=undetermined.
    #
    # Re-running is safe and idempotent: a hit already marked superseded stays
    # superseded, and a pair that still does not overlap is still left alone.
    for cluster in clusters:
        for family in families:
            own = [h for h in cluster.hits if h.family_key == family.key]
            if len(own) < 2:
                continue
            resolved_own, own_events = resolve_idiomorph_overlaps(
                own, family, min_overlap_fraction=idiomorph_overlap_fraction
            )
            # A pair the first pass already resolved comes back unchanged.
            # Recording it again reported every such event twice.
            already = idiomorph_events_by_cluster.get(id(cluster), [])
            own_events = [e for e in own_events if e not in already]
            if not own_events:
                continue
            by_id = dict(zip((id(h) for h in own), resolved_own))
            cluster.hits[:] = [by_id.get(id(h), h) for h in cluster.hits]
            idiomorph_events_by_cluster.setdefault(id(cluster), []).extend(own_events)
            resolutions_with_family.extend((family, e) for e in own_events)
            if evidence_diagnostics_path is not None:
                for event in own_events:
                    _write_idiomorph_diagnostics(
                        evidence_diagnostics_path, event, family,
                        run_id=run_id, genome_id=genome_id,
                    )

    # Cross-contig fragment merging is OFF by default. Curator's ruling,
    # 2026-09-20, after Basidiomycota order testing found a real Cryptococcus
    # deneoformans MAT locus -- both genes matching at 100% identity --
    # reported as one merged `fragmented` call whose top-level contig/start/end
    # named a chromosome holding almost none of its own evidence (a 3-segment
    # merge over 2 chromosomes, one span 849 kb), giving idiomorph:undetermined
    # and confidence: low for a strain whose mating type was not in doubt.
    #
    # Disabled, a genuinely fragmented locus is not lost: each contig's own
    # partial evidence is still reported by the normal per-cluster loop below,
    # just as separate, honest, correctly-coordinated calls rather than one
    # call that can misname itself. The merge machinery itself is unchanged
    # and still exercised by its own tests; only the default flipped.
    fragmented_segments: dict[FamilyKey, list[GeneCluster]] = {}
    if allow_cross_contig_fragments:
        for family in families:
            segments = _fragmented_family_segments(clusters, family)
            if segments:
                fragmented_segments[family.key] = segments

    def _any_gene_unpolished(cluster_ids: set[int], family_key: FamilyKey) -> bool:
        """True when ANY of this family's own genes, in these exact clusters,
        was localized but neither polishing tool could model it.

        Scoped by (cluster, family, gene): a gene left unpolished in one cluster
        must never cap the tier of an independent cluster of the same family,
        nor of another family that merely shares a gene name. Agreement
        (`polished_agree` vs `polished_disagree`) is deliberately NOT consulted
        -- both are "polished" as far as tiering is concerned.
        """
        return any(
            outcome.status == STATUS_UNPOLISHED
            for (cluster_id, key, _gene_name), outcome in polish_by.items()
            if key == family_key and cluster_id in cluster_ids
        )

    def _modelled_gene_count(
        member_clusters: list[GeneCluster], family_key: FamilyKey
    ) -> int:
        return len(_modelled_gene_names(member_clusters, family_key))

    def _modelled_gene_names(
        member_clusters: list[GeneCluster], family_key: FamilyKey
    ) -> frozenset[str]:
        """How many DISTINCT genes of this family, in these exact clusters, rest
        on a real gene model rather than on a bare alignment.

        A gene counts when either:

        * a polishing tool modelled it (`polished_agree`, `polished_disagree`
          or `polished_single`), or
        * it came from the annotated fast path, i.e. it has a
          `diamond_proteome` hit here -- a gene model somebody already called,
          which was never put to the tools precisely BECAUSE it needs no
          refining.

        Counting the fast path is not a loosening, it is the whole reason this
        is not called `polished_gene_count`. Gating on polish alone would score
        every gene of a fully annotated genome as unmodelled and withhold its
        real MAT locus -- the annotated ZygoLife and BFD-proteome workflows
        would report nothing at all. The measurement behind the bar came from
        genome-only runs, where no diamond hit exists and the two definitions
        coincide, so it does not speak to that case either way.

        What never counts is `STATUS_UNPOLISHED`: a localized HSP that BOTH
        tools were given a window for and neither could turn into a gene. That
        is the population the bar exists to remove.

        Nor does a gene whose every hit here is superseded. A RESCUED model is
        recorded as modelled before the post-polish idiomorph pass runs, and
        that pass can then find it is the losing half of a sexM/sexP-style
        pair -- the same physical gene as the winner. Counting it let one gene
        clear a bar of two (25 of 3,101 Saccharomyces loci were over-counted).

        Counted per distinct gene NAME, not per hit: a cluster routinely holds
        many HSPs of one gene, and three fragments of one alpha-box must not
        add up to the bar on their own. Scoped by (cluster, family) exactly as
        `_any_gene_unpolished` is, and for the same reasons.
        """
        cluster_ids = {id(c) for c in member_clusters}
        live = [
            h for c in member_clusters for h in c.hits
            if h.family_key == family_key and h.superseded_by is None
        ]
        modelled = {
            gene_name
            for (cluster_id, key, gene_name), outcome in polish_by.items()
            if key == family_key and cluster_id in cluster_ids
            and outcome.status in (STATUS_AGREE, STATUS_DISAGREE, STATUS_SINGLE)
        } & {h.gene_name for h in live}
        modelled |= {h.gene_name for h in live if h.method == "diamond_proteome"}
        return frozenset(modelled)

    def _build(
        score: FamilyScore,
        member_clusters: list[GeneCluster],
        scores_in_context: list[FamilyScore],
        fragmented: bool,
    ) -> DetectionResult:
        family = families_by_key[score.family_key]
        ambiguous = is_ambiguous(scores_in_context, floor=ambiguity_floor)
        tier = assign_tier(
            score, family, member_clusters[0],
            any_gene_unpolished=_any_gene_unpolished(
                {id(c) for c in member_clusters}, score.family_key
            ),
            fragmented=fragmented,
        )
        # The idiomorph calls behind this result, and the narrowest of them.
        # A thin margin does NOT withhold the call -- refusing would cost real
        # detections, and the rule is 23/23 correct on the ground-truth set
        # even at a 2.34-point separation. It caps the tier instead, so a
        # close call is visible at a glance rather than only to whoever opens
        # the diagnostics.
        own_resolutions = [
            event
            for cluster in member_clusters
            for event in idiomorph_events_by_cluster.get(id(cluster), ())
        ]
        # The margin of the call that was ACTUALLY MADE. Since the idiomorph is
        # decided by the bitscore vote rather than by overlap resolution, the
        # margin must come from the same ranking -- reporting the resolution's
        # identity margin here described a different mechanism and a different
        # quantity. Measured on Actinomucor sp. NRRL A-23671: this field read
        # 1.077 (identity) while the decision turned on a bitscore gap of 4.6.
        # Curator's ruling 2026-09-21: report the margin.
        #
        # Falls back to the resolution margin when the vote produced no ranking
        # (no idiomorph-restricted gene found), so nothing that used to be
        # reported is lost. Each resolution's own margin is unchanged and still
        # carried per event in `idiomorph_resolutions`.
        own_vote = idiomorph_candidates(
            family, score.genes_found,
            [h for c in member_clusters for h in _own_live_hits(c, family.key)],
        )
        idiomorph_margin = idiomorph_margin_from_vote(own_vote)
        if idiomorph_margin is None and own_resolutions:
            idiomorph_margin = min(e.margin for e in own_resolutions)
        if (
            idiomorph_margin is not None
            and idiomorph_margin < family.min_idiomorph_margin
            and tier == "high"
        ):
            tier = "medium"
        short_genes = short_orf_by_family.get(score.family_key, set())
        evidence = _gene_evidence(member_clusters, score.family_key, polish_by)
        # Segments are widened to cover this result's own gene evidence, so a
        # polished or rescued gene can never fall outside the locus segment
        # that reports it.
        segments = _segments_for(member_clusters, contig_lengths, evidence)
        # Classified BEFORE the result is built, because the class decides the
        # confidence cap: a `partial_locus` may never be `high`.
        partial_strength = is_partial_strength(
            score.fraction_found, ambiguity_floor, relaxed=False
        )
        locus_class = apply_partial_locus(
            classify_locus(member_clusters[0], family, full_length_models=full_length_models(
                polish_by, {id(c) for c in member_clusters}, score.family_key,
                record_protein_lengths,
                {g for c in member_clusters
                 for e in idiomorph_events_by_cluster.get(id(c), ())
                 for g in (e.winner, e.loser)},
            )),
            fraction_found=score.fraction_found,
            ambiguity_floor=ambiguity_floor,
            relaxed=False,
        )
        polished_genes = _modelled_gene_count(member_clusters, score.family_key)
        if polished_genes == 0:
            # Curator's rulings, 2026-09-22. (1) "without polishing it is low":
            # a locus whose every gene is a raw tblastn HSP that no tool could
            # model is not a medium-confidence call. Measured: 39,838 such loci
            # across the Pezizomycotina panels, mean identity 32.8%, and not
            # one of them is high-confidence today -- so nothing is demoted
            # from high by this.
            tier = "low"
            # (2) "same as 1, must have polished": `mat_locus` is the strongest
            # claim this pipeline makes, and 3,388 loci were making it on
            # entirely unpolished evidence. Demote to `partial_locus`, which is
            # already the strength class `apply_partial_locus` uses, rather
            # than inventing another. The composition classes
            # (`idiomorph_gene_only`, `flanking_gene_only`,
            # `homothallic_candidate`) are statements about WHICH genes are
            # present and are left alone, exactly as the floor-tie rule leaves
            # them.
            if locus_class == LOCUS_CLASS_MAT:
                locus_class = LOCUS_CLASS_PARTIAL
        return DetectionResult(
            polished_genes=polished_genes,
            family_key=score.family_key,
            contig=segments[0].contig,
            start=segments[0].start,
            # A fragmented call spans more than one contig, so there is no valid
            # single (contig, start, end); the top-level coordinates name the
            # first segment and `segments` carries the rest.
            end=segments[0].end,
            # Capped on STRENGTH, not on the resulting class: a
            # `homothallic_candidate` that only tied the floor keeps its class
            # (it is a statement about which genes are present) but must not
            # keep `high` (that is a statement about how sure we are).
            confidence=cap_at_medium(tier) if partial_strength else tier,
            idiomorph=assign_idiomorph(
                family, score.genes_found,
                [h for c in member_clusters for h in _own_live_hits(c, family.key)],
            ),
            idiomorph_candidates=own_vote,
            ambiguous_with=(
                [
                    s.family_key
                    for s in scores_in_context
                    if s.family_key != score.family_key and s.fraction_found >= ambiguity_floor
                ]
                if ambiguous
                else []
            ),
            genes_found=score.genes_found,
            genes_missing=[g for g in score.genes_missing if g not in short_genes],
            fragmented=fragmented,
            # Two reasons a gene could not be searched for, reported as one
            # list because they mean the same thing to a reader: scoring
            # already excluded the genes with no reference protein at all,
            # and `short_genes` names those whose only reference is too short
            # to localize reliably.
            genes_not_searchable=sorted(
                set(score.genes_not_searchable)
                | {g for g in score.genes_missing if g in short_genes}
            ),
            segments=segments,
            gene_evidence=evidence,
            reference_records=sorted({e.reference_record_id for e in evidence}),
            idiomorph_margin=idiomorph_margin,
            idiomorph_resolutions=own_resolutions,
            locus_class=locus_class,
            # Measured across the whole call, including every segment of a
            # fragmented one -- the span a reader sees in the report is the
            # thing being judged.
            span_exceeds_plausible_bound=span_exceeds_plausible_bound(
                min(s.start for s in segments), max(s.end for s in segments),
                family,
            ),
        )

    results: list[DetectionResult] = []
    # best (fraction, score, clusters) seen per family, for "not detected" reporting
    best_attempt: dict[FamilyKey, FamilyScore] = {}
    # Cluster ids actually reported as segments of a fragmented multi-segment
    # call, keyed per family. Only these specific clusters are suppressed from
    # the per-cluster loop below -- a family judged fragmented may still have
    # a separate, independent, above-floor cluster elsewhere in the genome
    # (e.g. a genuine second locus from gene duplication), and that cluster
    # must still go through normal per-cluster reporting rather than being
    # dropped just because it shares a family_key with the fragmented call.
    fragmented_reported_cluster_ids: dict[FamilyKey, set[int]] = {}

    for family_key, member_clusters in fragmented_segments.items():
        merged = GeneCluster(
            member_clusters[0].contig,
            member_clusters[0].start,
            member_clusters[0].end,
            [h for c in member_clusters for h in c.hits],
        )
        scores = score_cluster(merged, families, searchable_genes=searchable_genes)
        score = next((s for s in scores if s.family_key == family_key), None)
        if score is None:
            continue
        best_attempt[family_key] = score
        if score.fraction_found < ambiguity_floor and not is_ambiguous(scores, floor=ambiguity_floor):
            continue
        results.append(_build(score, member_clusters, scores, fragmented=True))
        fragmented_reported_cluster_ids[family_key] = {id(c) for c in member_clusters}

    for cluster in clusters:
        scores = score_cluster(cluster, families, searchable_genes=searchable_genes)
        ambiguous = is_ambiguous(scores, floor=ambiguity_floor)
        for score in scores:
            if id(cluster) in fragmented_reported_cluster_ids.get(score.family_key, ()):
                continue  # this exact cluster was already reported as a segment of the fragmented call
            previous = best_attempt.get(score.family_key)
            if previous is None or score.fraction_found > previous.fraction_found:
                best_attempt[score.family_key] = score
            if score.fraction_found < ambiguity_floor and not ambiguous:
                continue
            results.append(_build(score, [cluster], scores, fragmented=False))

    # The relaxed second pass, and ONLY when the strict pass found nothing
    # anywhere in this genome. Gating on the whole genome rather than per
    # family is the curator's ruling and the conservative reading: a genome
    # with any confident call is left exactly as it was, so the relaxed bar
    # can never dilute a run that already worked.
    if not results and relaxed_second_pass:
        results = _relaxed_results(
            clusters, families,
            searchable_genes=searchable_genes,
            evidence_floor=evidence_floor,
            ambiguity_floor=ambiguity_floor,
            polish_by=polish_by,
            contig_lengths=contig_lengths,
        )
        if results:
            logger.info(
                "strict pass found no locus; %d admitted by the relaxed pass "
                "(>=2 distinct genes incl. a core gene), capped at medium",
                len(results),
            )

    # Curator's ruling (3), 2026-09-22: "yes require more than 1". A cluster
    # must carry at least MIN_POLISHED_GENES polished gene models to be
    # REPORTED. Applied here, after the relaxed pass, so the relaxed trigger
    # ("only when the strict pass found nothing") still sees the unfiltered
    # strict result and its behaviour is unchanged -- and so the bar applies
    # to relaxed calls too, which is the point of having it.
    #
    # Nothing is destroyed: the per-candidate rows are already on disk in the
    # evidence-diagnostics stream, a family whose only calls were withheld
    # still gets a `not_detected` entry naming the reason, and the count is
    # carried on the outcome.
    suppressed = [r for r in results if r.polished_genes < min_polished_genes]
    results = [r for r in results if r.polished_genes >= min_polished_genes]
    if suppressed:
        logger.info(
            "withheld %d locus/loci carrying fewer than %d modelled genes",
            len(suppressed), min_polished_genes,
        )

    reported = {r.family_key for r in results}
    #: family -> the best withheld call for it, so `not_detected` can say that
    #: something WAS built and why it did not survive, rather than falling
    #: through to the generic below-the-floor wording.
    suppressed_best: dict[FamilyKey, DetectionResult] = {}
    for r in suppressed:
        best = suppressed_best.get(r.family_key)
        if best is None or r.polished_genes > best.polished_genes:
            suppressed_best[r.family_key] = r
    not_detected: list[NotDetectedFamily] = []
    for family in families:
        if family.key in reported:
            continue
        short_genes = short_orf_by_family.get(family.key, set())
        score = best_attempt.get(family.key)
        withheld = suppressed_best.get(family.key)
        if withheld is not None:
            modelled = (
                "none of its genes could be modelled" if withheld.polished_genes == 0
                else f"only {withheld.polished_genes} of its genes could be modelled"
            )
            not_detected.append(NotDetectedFamily(
                family_key=family.key,
                reason=(
                    f"best cluster carried {withheld.polished_genes} modelled "
                    f"gene(s), below the {min_polished_genes} required to report a "
                    f"locus; {modelled}"
                ),
                best_fraction_found=score.fraction_found if score else 0.0,
                genes_found=list(withheld.genes_found),
                genes_missing=[g for g in withheld.genes_missing if g not in short_genes],
                genes_not_searchable=sorted(short_genes),
            ))
            continue
        if score is None and reference_fasta.exists() and not searchable_genes.get(family.key):
            # Nothing to search with: every record of the family was withheld
            # (a holdout) or none is curated. "No hits" would blame the genome.
            # A MISSING file is "no information" (see searchable_genes_by_family),
            # not an empty set, so it keeps the generic wording below.
            not_detected.append(NotDetectedFamily(
                family_key=family.key,
                reason="no reference protein for this family is in the search set "
                       "(none curated, or all withheld)",
                best_fraction_found=0.0,
                genes_missing=[g["name"] for g in family.genes if g["name"] not in short_genes],
                genes_not_searchable=sorted(short_genes),
            ))
        elif score is None:
            not_detected.append(NotDetectedFamily(
                family_key=family.key,
                reason="no reference-protein hits found for this family in this genome",
                best_fraction_found=0.0,
                genes_missing=[g["name"] for g in family.genes if g["name"] not in short_genes],
                genes_not_searchable=sorted(short_genes),
            ))
        else:
            not_detected.append(NotDetectedFamily(
                family_key=family.key,
                reason=(
                    f"best cluster matched {score.fraction_found:.2f} of this family's expected "
                    f"genes, below the ambiguity floor of {ambiguity_floor:.2f}"
                ),
                best_fraction_found=score.fraction_found,
                genes_found=score.genes_found,
                genes_missing=[g for g in score.genes_missing if g not in short_genes],
                genes_not_searchable=[g for g in score.genes_missing if g in short_genes],
            ))

    return DetectionOutcome(
        results=results,
        not_detected=not_detected,
        families_attempted=[f.key for f in families],
        routing_mode=routing.routing_mode,
        routing_error=routing.routing_error,
        genetic_code=genetic_code,
        genetic_code_error=genetic_code_error,
        suppressed_unpolished=len(suppressed),
        suppressed_loci=suppressed,
    )
