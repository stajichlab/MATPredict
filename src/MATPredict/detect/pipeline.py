"""Orchestrates routing -> search -> clustering -> scoring -> tiering -> idiomorph assignment.

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
from MATPredict.detect.scoring import FamilyScore, is_ambiguous, score_cluster
from MATPredict.detect.search import SearchHit, search_fast_path, search_genomic
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


def _gene_evidence(hits: list[SearchHit], family_key: FamilyKey) -> list[GeneEvidence]:
    """Best (highest-identity) hit per gene name for one family, as report-ready evidence."""
    best: dict[str, SearchHit] = {}
    for hit in hits:
        if hit.family_key != family_key:
            continue
        current = best.get(hit.gene_name)
        if current is None or hit.identity > current.identity:
            best[hit.gene_name] = hit
    return [
        GeneEvidence(
            gene_name=h.gene_name, role=h.role, contig=h.contig, start=h.start, end=h.end,
            strand=h.strand, identity=h.identity, coverage=h.coverage,
            reference_record_id=h.reference_record_id, method=h.method,
        )
        for h in sorted(best.values(), key=lambda h: (h.contig, h.start))
    ]


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
    search_genomic: Callable = search_genomic,
    max_gap: int = 25_000,
    ambiguity_floor: float = 0.5,
    short_orf_aa_floor: int = 60,
) -> DetectionOutcome:
    families = route(taxid, load_all_families(db_root))
    record_families = load_record_families(db_root)
    short_orf_by_family = _short_orf_genes(
        db_root, families, record_families, short_orf_aa_floor
    )
    families_by_key = {f.key: f for f in families}

    hits: list[SearchHit] = []
    if proteome_fasta is not None:
        hits.extend(search_fast_path(proteome_fasta, families, reference_fasta, record_families))
    else:
        hits.extend(search_genomic(genome_fasta, families, reference_fasta, record_families))

    # Unconditional second pass, part 1 -- whole-genome, for families with NO
    # hits at all. The spec requires the genomic second pass to fire for "any
    # expected core_MAT gene not found among the proteome hits ... unconditional,
    # not gated on whether flanking genes were found first". A family whose genes
    # are ALL missing from a supplied proteome annotation (the mfa1-style blind
    # spot this pipeline exists for) has no cluster to anchor a window on, so
    # anchoring the second pass to an existing cluster would skip exactly the
    # case that matters most. Those families get a relaxed search against the
    # whole genome instead.
    #
    # All zero-hit families are searched in ONE batched exonerate call rather
    # than one call each: a whole-genome protein2genome run is expensive, and
    # the exhaustive (no-taxid) path routes every family in the database.
    families_with_hits = {h.family_key for h in hits}
    missing_families = [f for f in families if f.key not in families_with_hits]
    genome_wide_second_pass: set[FamilyKey] = set()
    if missing_families:
        missing_keys = {f.key for f in missing_families}
        rescued = [
            h
            for h in search_genomic(
                genome_fasta, missing_families, reference_fasta, record_families, relaxed=True
            )
            if h.family_key in missing_keys
        ]
        genome_wide_second_pass = {h.family_key for h in rescued}
        hits.extend(rescued)

    clusters = cluster_hits(hits, max_gap=max_gap)

    # Unconditional second pass, part 2 -- windowed, for any family that has a
    # foothold in a cluster but is still missing one of its own core_MAT genes.
    # Each family's missing genes are checked strictly against that same
    # family's own hits (see _missing_core_genes), so one family's presence
    # never masks or substitutes for another family's absence.
    #
    # second_pass_used_for is keyed by (id(cluster), family.key), NOT by
    # family.key alone. The same family can have multiple independent
    # spatial clusters in one genome (e.g. gene-duplication / multi-allele
    # co-occurrence at MAT loci), and whether the relaxed second pass was
    # needed in one cluster must never leak into the tiering of an unrelated
    # cluster for the same family.
    second_pass_used_for: set[tuple[int, FamilyKey]] = set()
    for cluster in clusters:
        for family in _families_with_a_foothold(cluster, families):
            if not _missing_core_genes(cluster, family):
                continue
            relaxed_hits = search_genomic(
                genome_fasta, [family], reference_fasta, record_families, relaxed=True,
                window=(cluster.contig, cluster.start, cluster.end),
            )
            # Only accept hits that actually belong to this specific family --
            # a defensive filter in case a stubbed/real search_genomic ever
            # returns hits keyed to a different family_key than requested.
            own_hits = [h for h in relaxed_hits if h.family_key == family.key]
            if own_hits:
                second_pass_used_for.add((id(cluster), family.key))
                cluster.hits.extend(own_hits)

    fragmented_segments: dict[FamilyKey, list[GeneCluster]] = {}
    for family in families:
        segments = _fragmented_family_segments(clusters, family)
        if segments:
            fragmented_segments[family.key] = segments

    # Read contig lengths unconditionally so contig_edge_distance is populated
    # (or left None on an unreadable FASTA) consistently for every segment of
    # every result, independent of whether some OTHER family in this same run
    # happened to be fragmented. Gating this on `fragmented_segments` used to
    # make two otherwise-identical single-contig runs disagree on
    # contig_edge_distance purely because of an unrelated family elsewhere in
    # the genome.
    contig_lengths = _contig_lengths(genome_fasta)

    def _second_pass_used(cluster_ids: list[int], family_key: FamilyKey) -> bool:
        return family_key in genome_wide_second_pass or any(
            (cid, family_key) in second_pass_used_for for cid in cluster_ids
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
            second_pass_used=_second_pass_used([id(c) for c in member_clusters], score.family_key),
            fragmented=fragmented,
        )
        short_genes = short_orf_by_family.get(score.family_key, set())
        all_hits = [h for c in member_clusters for h in c.hits]
        evidence = _gene_evidence(all_hits, score.family_key)
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
                reason="no reference-protein hits found for this family in this genome, "
                       "including after the relaxed whole-genome second pass",
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
