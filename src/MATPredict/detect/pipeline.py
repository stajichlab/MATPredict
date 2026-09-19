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

Fragmented assemblies (spec section 6) are handled to the extent described in
`_fragmented_family_segments`: a family whose expected `core_MAT` genes are
split across clusters on *different contigs*, with no single cluster carrying
them all, is reported once as a multi-segment call with `fragmented=True`
(which downgrades its confidence tier by one, per `tiering.assign_tier`).
`contig_edge_distance` is populated per segment when the genome FASTA can be
read; it is left as None otherwise. What is deliberately NOT attempted here is
reconstructing locus order or the intervening sequence across segments -- a
multi-segment call reports the segments it found, nothing more.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from MATPredict.detect.clustering import GeneCluster, cluster_hits
from MATPredict.detect.family_registry import (
    Family,
    FamilyKey,
    load_all_families,
    load_record_families,
    route,
)
from MATPredict.detect.idiomorph import assign_idiomorph
from MATPredict.detect.polish import (
    STATUS_DISAGREE,
    STATUS_NOT_POLISH_CANDIDATE,
    STATUS_UNPOLISHED,
    PolishModel,
    PolishOutcome,
    classify,
)
from MATPredict.detect.scoring import FamilyScore, is_ambiguous, score_cluster
from MATPredict.detect.search import (
    SearchHit,
    polish_with_exonerate,
    polish_with_miniprot,
    search_fast_path,
    search_localize,
)
from MATPredict.detect.tiering import assign_tier


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


def _missing_core_genes(cluster: GeneCluster, family: Family) -> set[str]:
    """A single family's own core_MAT genes not yet found among that same
    family's own hits in this cluster -- never checked against another
    family's hits or expected genes (see Task 6's cross-family pooling bug)."""
    core_genes = {g["name"] for g in family.genes if g["role"] == "core_MAT"}
    found_genes = {h.gene_name for h in cluster.hits if h.family_key == family.key}
    return core_genes - found_genes


def _curated_protein_lengths(
    db_root: Path,
    families: list[Family],
    record_families: dict[str, FamilyKey],
) -> dict[tuple[FamilyKey, str], int]:
    """`(family_key, gene_name)` -> the LONGEST curated reference protein, in aa.

    Keyed per `(phylum, locus_name)` family, never by bare gene name: gene names
    are reused across families (`sla2` is both `Ascomycota:MATsc`'s and
    `Ascomycota:MATyl`'s), so one family's entry must never answer for another's.
    See `_short_orf_genes` for why the MAXIMUM is the right summary.

    Reads db/**/proteins.faa (matching gff_export.write_proteins_fasta's header
    form `>{record_id}|gene_index={n}|name={name}|role={role}`, the same raw
    per-record files reference_fasta.py concatenates -- not its rewritten output).
    """
    expected_by_family = {f.key: {g["name"] for g in f.genes} for f in families}
    longest: dict[tuple[FamilyKey, str], int] = {}

    for faa in db_root.glob("*/*/*/proteins.faa"):
        if faa.relative_to(db_root).parts[0] == "candidates":
            continue
        for chunk in faa.read_text().split(">")[1:]:
            header, _, seq = chunk.partition("\n")
            head_fields = header.split("|")
            record_id = head_fields[0]
            family_key = record_families.get(record_id)
            if family_key is None or family_key not in expected_by_family:
                continue
            parts = dict(p.split("=", 1) for p in head_fields[1:] if "=" in p)
            name = parts.get("name")
            if name not in expected_by_family[family_key]:
                continue
            length = len(seq.strip().replace("\n", ""))
            key = (family_key, name)
            longest[key] = max(longest.get(key, 0), length)
    return longest


def _short_orf_genes(
    db_root: Path,
    families: list[Family],
    record_families: dict[str, FamilyKey],
    floor_aa: int,
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
    longest = _curated_protein_lengths(db_root, families, record_families)

    short_by_family: dict[FamilyKey, set[str]] = {}
    for (family_key, name), length in longest.items():
        if length < floor_aa:
            short_by_family.setdefault(family_key, set()).add(name)
    return short_by_family


@dataclass(frozen=True)
class EvidenceFloor:
    """Minimum evidence a family must clear in a cluster before it is admitted
    to the (expensive, per-gene, two-subprocess-per-tool) Stage 2 polish loop.

    Defaults exactly reproduce the pre-existing `_families_with_a_foothold`
    behavior (any single hit of any role, any identity) -- this gate changes
    nothing until a caller sets a stricter floor. See run_pipeline's own
    `evidence_floor`/`evidence_diagnostics_path` parameters and this task's
    plan notes for why no non-default value is set here yet: this rollout's
    own real data shows a naive identity or hit-count cutoff is not yet safely
    separable from real signal, so tuning these is left to a calibration pass
    once db/Ascomycota/order.yml's taxonomic_scope fix (a separate task) lets
    most genomes run with correctly-narrowed routing instead of the
    exhaustive fallback that produced this rollout's noisy example.
    """

    min_hits: int = 1
    min_identity: float | None = None
    require_core_role: bool = False


def _families_meeting_evidence_floor(
    cluster: GeneCluster, families: list[Family], floor: EvidenceFloor
) -> list[Family]:
    """Families whose OWN hits in this cluster clear `floor` -- generalizes
    the retired `_families_with_a_foothold` (which was exactly
    `_families_meeting_evidence_floor(cluster, families, EvidenceFloor())`)."""
    admitted = []
    for family in families:
        own_hits = [h for h in cluster.hits if h.family_key == family.key]
        if len(own_hits) < floor.min_hits:
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


def _write_evidence_diagnostics(
    out_path: Path, cluster: GeneCluster, family: Family, admitted: bool
) -> None:
    """Append one JSON line describing this (cluster, family) admission
    decision -- the real calibration dataset `EvidenceFloor`'s docstring
    refers to. Never raises on a write failure; diagnostics are best-effort
    and must never abort a real detection run."""
    own_hits = [h for h in cluster.hits if h.family_key == family.key]
    row = {
        "family": f"{family.key.phylum}:{family.key.locus_name}",
        "contig": cluster.contig, "cluster_start": cluster.start, "cluster_end": cluster.end,
        "gene_count": len({h.gene_name for h in own_hits}),
        "roles": sorted({h.role for h in own_hits}),
        "best_identity": max((h.identity for h in own_hits), default=None),
        "admitted": admitted,
    }
    try:
        with out_path.open("a") as f:
            f.write(json.dumps(row) + "\n")
    except OSError:
        pass


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
        for hit in cluster.hits:
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
    max_gap: int = 25_000,
    ambiguity_floor: float = 0.5,
    short_orf_aa_floor: int = 60,
    window_protein_length_multiple: float = DEFAULT_WINDOW_PROTEIN_LENGTH_MULTIPLE,
    window_max_intron_bp: int = DEFAULT_WINDOW_MAX_INTRON_BP,
    polish_tolerance_bp: int = 10,
    evidence_floor: EvidenceFloor = EvidenceFloor(),
    evidence_diagnostics_path: Path | None = None,
) -> DetectionOutcome:
    families = route(taxid, load_all_families(db_root))
    record_families = load_record_families(db_root)
    protein_lengths = _curated_protein_lengths(db_root, families, record_families)
    short_orf_by_family = _short_orf_genes(
        db_root, families, record_families, short_orf_aa_floor
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
                    genome_fasta, rescue_scope.families, reference_fasta, record_families
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
        localized = search_localize(genome_fasta, families, reference_fasta, record_families)
        localized_hit_ids.update(id(h) for h in localized)
        hits.extend(localized)

    clusters = cluster_hits(hits, max_gap=max_gap)

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
            for family in _families_meeting_evidence_floor(cluster, families, EvidenceFloor()):
                _write_evidence_diagnostics(
                    evidence_diagnostics_path, cluster, family,
                    admitted=family.key in admitted_keys,
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
                exonerate_model = polish_with_exonerate(
                    genome_fasta=genome_fasta, family=family, gene_name=gene_name,
                    reference_fasta=reference_fasta, record_families=record_families,
                    window=window,
                )
                miniprot_model = polish_with_miniprot(
                    genome_fasta=genome_fasta, family=family, gene_name=gene_name,
                    reference_fasta=reference_fasta, record_families=record_families,
                    window=window,
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

    fragmented_segments: dict[FamilyKey, list[GeneCluster]] = {}
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
        short_genes = short_orf_by_family.get(score.family_key, set())
        evidence = _gene_evidence(member_clusters, score.family_key, polish_by)
        # Segments are widened to cover this result's own gene evidence, so a
        # polished or rescued gene can never fall outside the locus segment
        # that reports it.
        segments = _segments_for(member_clusters, contig_lengths, evidence)
        return DetectionResult(
            family_key=score.family_key,
            contig=segments[0].contig,
            start=segments[0].start,
            # A fragmented call spans more than one contig, so there is no valid
            # single (contig, start, end); the top-level coordinates name the
            # first segment and `segments` carries the rest.
            end=segments[0].end,
            confidence=tier,
            idiomorph=assign_idiomorph(family, score.genes_found),
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
            genes_not_searchable=[g for g in score.genes_missing if g in short_genes],
            segments=segments,
            gene_evidence=evidence,
            reference_records=sorted({e.reference_record_id for e in evidence}),
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
        scores = score_cluster(merged, families)
        score = next((s for s in scores if s.family_key == family_key), None)
        if score is None:
            continue
        best_attempt[family_key] = score
        if score.fraction_found < ambiguity_floor and not is_ambiguous(scores, floor=ambiguity_floor):
            continue
        results.append(_build(score, member_clusters, scores, fragmented=True))
        fragmented_reported_cluster_ids[family_key] = {id(c) for c in member_clusters}

    for cluster in clusters:
        scores = score_cluster(cluster, families)
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

    reported = {r.family_key for r in results}
    not_detected: list[NotDetectedFamily] = []
    for family in families:
        if family.key in reported:
            continue
        short_genes = short_orf_by_family.get(family.key, set())
        score = best_attempt.get(family.key)
        if score is None:
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
    )
