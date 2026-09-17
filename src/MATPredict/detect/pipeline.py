"""Orchestrates routing -> search -> clustering -> scoring -> tiering -> idiomorph assignment."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

from MATPredict.detect.clustering import GeneCluster, cluster_hits
from MATPredict.detect.family_registry import Family, FamilyKey, load_all_families, route
from MATPredict.detect.idiomorph import assign_idiomorph
from MATPredict.detect.scoring import is_ambiguous, score_cluster
from MATPredict.detect.search import SearchHit, search_fast_path, search_genomic
from MATPredict.detect.tiering import assign_tier


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


def _missing_core_genes(cluster: GeneCluster, family: Family) -> set[str]:
    """A single family's own core_MAT genes not yet found among that same
    family's own hits in this cluster -- never checked against another
    family's hits or expected genes (see Task 6's cross-family pooling bug)."""
    core_genes = {g["name"] for g in family.genes if g["role"] == "core_MAT"}
    found_genes = {h.gene_name for h in cluster.hits if h.family_key == family.key}
    return core_genes - found_genes


def _short_orf_genes(db_root: Path, families: list[Family], floor_aa: int) -> set[str]:
    """Gene names whose curated reference protein length is below floor_aa,
    scanned once from db/**/proteins.faa (matching gff_export.write_proteins_fasta's
    header form `>{record_id}|gene_index={n}|name={name}|role={role}`, the same raw
    per-record files reference_fasta.py concatenates -- not its rewritten output)."""
    short_genes: set[str] = set()
    expected = {g["name"] for f in families for g in f.genes}
    for faa in db_root.glob("*/*/*/proteins.faa"):
        for chunk in faa.read_text().split(">")[1:]:
            header, _, seq = chunk.partition("\n")
            parts = dict(p.split("=", 1) for p in header.split("|")[1:] if "=" in p)
            name = parts.get("name")
            if name in expected and len(seq.strip()) < floor_aa:
                short_genes.add(name)
    return short_genes


def _families_with_a_foothold(cluster: GeneCluster, families: list[Family]) -> list[Family]:
    """Families that already have at least one hit of their own in this
    cluster -- the second pass only ever re-searches for a family's own
    missing core gene, keyed strictly by that family's family_key."""
    families_with_hits = {h.family_key for h in cluster.hits}
    return [f for f in families if f.key in families_with_hits]


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
) -> list[DetectionResult]:
    families = route(taxid, load_all_families(db_root))
    short_orf_genes = _short_orf_genes(db_root, families, short_orf_aa_floor)

    hits: list[SearchHit] = []
    if proteome_fasta is not None:
        hits.extend(search_fast_path(proteome_fasta, families, reference_fasta))
    else:
        hits.extend(search_genomic(genome_fasta, families, reference_fasta))

    clusters = cluster_hits(hits, max_gap=max_gap)

    # Unconditional genomic re-search for any core_MAT gene missing from any
    # family that already has a foothold in a cluster -- never gated on
    # whether flanking genes were found (spec section "Search"). Each
    # family's missing genes are checked strictly against that same
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
                genome_fasta, [family], reference_fasta, relaxed=True,
                window=(cluster.contig, cluster.start, cluster.end),
            )
            # Only accept hits that actually belong to this specific family --
            # a defensive filter in case a stubbed/real search_genomic ever
            # returns hits keyed to a different family_key than requested.
            own_hits = [h for h in relaxed_hits if h.family_key == family.key]
            if own_hits:
                second_pass_used_for.add((id(cluster), family.key))
                cluster.hits.extend(own_hits)

    results: list[DetectionResult] = []
    families_by_key = {f.key: f for f in families}
    for cluster in clusters:
        scores = score_cluster(cluster, families)
        if not scores:
            continue
        ambiguous = is_ambiguous(scores, floor=ambiguity_floor)
        for score in scores:
            if score.fraction_found < ambiguity_floor and not ambiguous:
                continue
            family = families_by_key[score.family_key]
            tier = assign_tier(
                score, family, cluster,
                second_pass_used=(id(cluster), score.family_key) in second_pass_used_for,
                fragmented=False,
            )
            idiomorph = assign_idiomorph(family, score.genes_found)
            ambiguous_with = (
                [s.family_key for s in scores if s.family_key != score.family_key and s.fraction_found >= ambiguity_floor]
                if ambiguous else []
            )
            genes_not_searchable = [g for g in score.genes_missing if g in short_orf_genes]
            genes_missing = [g for g in score.genes_missing if g not in short_orf_genes]
            results.append(DetectionResult(
                family_key=score.family_key, contig=cluster.contig, start=cluster.start, end=cluster.end,
                confidence=tier, idiomorph=idiomorph, ambiguous_with=ambiguous_with,
                genes_found=score.genes_found, genes_missing=genes_missing, fragmented=False,
                genes_not_searchable=genes_not_searchable,
            ))
    return results
