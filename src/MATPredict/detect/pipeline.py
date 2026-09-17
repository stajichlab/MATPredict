"""Orchestrates routing -> search -> clustering -> polishing -> scoring -> tiering -> idiomorph assignment.

Search is a two-stage *localize-then-polish* flow (see
docs/superpowers/specs/2026-09-17-mat-detection-search-localization-design.md):

* **Genome-only path** (no predicted proteome): one batched genome-wide
  `tblastn` localization (`search_localize`) over every routed family's
  reference proteins, clustered once by `cluster_hits`. Every gene that
  `tblastn` localized in a cluster is then polished.
* **Fast-path rescue** (a proteome was supplied): `search_fast_path` runs as
  before; a cluster/family that has a foothold but is still missing one of its
  own `core_MAT` genes skips localization entirely -- a rough location already
  exists -- and polishes that one gene directly against a window padded around
  the *existing* cluster's span.

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


def _families_with_a_foothold(cluster: GeneCluster, families: list[Family]) -> list[Family]:
    """Families that already have at least one hit of their own in this
    cluster -- the windowed second pass only ever re-searches for a family's own
    missing core gene, keyed strictly by that family's family_key."""
    families_with_hits = {h.family_key for h in cluster.hits}
    return [f for f in families if f.key in families_with_hits]


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
    model: PolishModel | None, family_key: FamilyKey, gene_name: str
) -> PolishModel | None:
    """A polished model, or None when it is not this family's own named gene.

    Both polishing wrappers already drop a model whose gene name or curated
    record does not match what was asked for, so this is defence in depth --
    kept because the retired windowed second pass carried the same filter, and
    because a model silently attributed to the wrong family is precisely the
    cross-family bug class this file has shipped before. Dropping it degrades
    the gene to `unpolished` (its raw localization hit still stands) rather than
    writing another family's coordinates into this family's evidence.
    """
    if model is None:
        return None
    if model.family_key != family_key or model.gene_name != gene_name:
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


def _gene_evidence(
    member_clusters: list[GeneCluster],
    family_key: FamilyKey,
    polish_by: dict[tuple[int, FamilyKey, str], PolishOutcome],
) -> list[GeneEvidence]:
    """Report-ready per-gene evidence for ONE family across its own clusters.

    Per gene, the canonical `PolishOutcome.canonical` model wins when that gene
    was polished in that cluster; a `STATUS_UNPOLISHED` gene (canonical is None)
    and any gene that was never a polish candidate at all -- e.g. a gene the
    diamond fast path already found -- fall back to that gene's best
    (highest-identity) raw `SearchHit`. That fallback is the only path by which
    a `GeneEvidence` is built from a `SearchHit` rather than a `PolishModel`.

    `polish_by` is keyed `(id(cluster), family_key, gene_name)` and is read here
    ONLY for this family's own genes in these exact clusters, so a polish result
    from another family, or from an unrelated cluster of this same family, can
    never be attributed to this call.
    """
    best: dict[str, GeneEvidence] = {}
    for cluster in member_clusters:
        raw_by_gene: dict[str, SearchHit] = {}
        for hit in cluster.hits:
            if hit.family_key != family_key:
                continue
            current = raw_by_gene.get(hit.gene_name)
            if current is None or hit.identity > current.identity:
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
                evidence = GeneEvidence(
                    gene_name=model.gene_name, role=model.role, contig=model.contig,
                    start=model.start, end=model.end, strand=model.strand,
                    identity=model.identity, coverage=None,
                    reference_record_id=model.reference_record_id, method=model.method,
                )
            else:
                hit = raw_by_gene.get(gene_name)
                if hit is None:
                    continue  # unpolished with nothing localized: no evidence to report
                evidence = GeneEvidence(
                    gene_name=hit.gene_name, role=hit.role, contig=hit.contig,
                    start=hit.start, end=hit.end, strand=hit.strand,
                    identity=hit.identity, coverage=hit.coverage,
                    reference_record_id=hit.reference_record_id, method=hit.method,
                )
            previous = best.get(gene_name)
            if previous is None or evidence.identity > previous.identity:
                best[gene_name] = evidence
    return sorted(best.values(), key=lambda e: (e.contig, e.start))


def _segments_for(
    clusters: list[GeneCluster], contig_lengths: dict[str, int]
) -> list[LocusSegment]:
    segments = []
    for cluster in sorted(clusters, key=lambda c: (c.contig, c.start)):
        length = contig_lengths.get(cluster.contig)
        edge = min(cluster.start - 1, length - cluster.end) if length else None
        segments.append(LocusSegment(cluster.contig, cluster.start, cluster.end, edge))
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
) -> DetectionOutcome:
    families = route(taxid, load_all_families(db_root))
    record_families = load_record_families(db_root)
    protein_lengths = _curated_protein_lengths(db_root, families, record_families)
    short_orf_by_family = _short_orf_genes(
        db_root, families, record_families, short_orf_aa_floor
    )
    families_by_key = {f.key: f for f in families}

    # Stage 0/1 -- one search per run. The genome-only path localizes with a
    # single batched genome-wide tblastn call covering every routed family; the
    # fast path uses the supplied proteome and does NOT localize at all.
    hits: list[SearchHit] = []
    if proteome_fasta is not None:
        hits.extend(search_fast_path(proteome_fasta, families, reference_fasta, record_families))
    else:
        hits.extend(search_localize(genome_fasta, families, reference_fasta, record_families))

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
        for family in _families_with_a_foothold(cluster, families):
            if proteome_fasta is None:
                # Genome-only path: every gene tblastn localized in THIS cluster
                # for THIS family needs its approximate HSP coordinates refined.
                genes_to_polish = sorted(
                    {h.gene_name for h in cluster.hits if h.family_key == family.key}
                )
                rescue = False
            else:
                # Fast-path rescue: only this family's own still-missing core
                # genes, checked strictly against this family's own hits in this
                # cluster (see _missing_core_genes), so one family's presence
                # never masks another family's absence.
                genes_to_polish = sorted(_missing_core_genes(cluster, family))
                rescue = True

            for gene_name in genes_to_polish:
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
                    _own_model(exonerate_model, family.key, gene_name),
                    _own_model(miniprot_model, family.key, gene_name),
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
        segments = _segments_for(member_clusters, contig_lengths)
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
