"""A MAT locus split across contigs by the assembly.

Curator's ruling 2026-09-27 (results/2026-09-27_rarrhizus_uncalled/NOTE.md):
a single MODELLED core gene at >= 95% identity lying near a contig end, whose
family's flanking genes are found on OTHER contigs at similar identity, is
reported as `partial_locus`, confidence `low`, flagged `split_locus`, instead
of being withheld by the two-modelled-gene bar.

Why: clusters never cross contigs (by design), so when an assembly breaks
between the core gene and its flanks, a real locus is left as one modelled
gene on one contig and cannot reach the bar of two. Measured case: 12
Rhizopus arrhizus/delemar GL-series assemblies carry sexP at 98.4% alone on a
contig 1-199 bp from its end (median 167), tptA/btbA on a second contig and
rnhA on a third -- all at 98.6-100%. Same failure mode as the split
Cryptococcus loci and the C. auris gap cases.

The defaults and why:

* `min_core_identity` 95 -- the ruling. A near-identical match to a curated
  record is what makes one gene strong enough to stand for the locus; the GL
  cases are 98.4%, HMG paralogs in the same genomes are ~40%.
* `max_edge_bp` 500 -- the GL cores end 1-199 bp from the contig end; 500 is
  2.5x the largest observed, so a slightly different break still qualifies,
  while a gene deep inside a contig (where the flanks would be expected on the
  same contig) does not.
* `min_flank_genes` 2 and `flank_identity_gap` 10 -- at least two distinct
  roster flanking genes, each at >= core identity - 10 points. One flank alone
  is too easy to meet by chance; the GL flanks sit within 1.6 points of the
  core. At least one of them must lie on a DIFFERENT contig from the core --
  otherwise the locus is simply unsplit and the ordinary rules apply.
* Paralog guards: the core hit must be the genome's best (highest bitscore)
  hit for its gene; each flank is represented by the genome's best hit for
  that flank gene; a flank inside a locus already reported (any family) does
  not count.

The idiomorph comes from `_build` as for any call, i.e. from the family's
HMM classifier on the single modelled core protein when it has one.
"""
from __future__ import annotations

from dataclasses import dataclass, fields

#: polish.STATUS_AGREE / STATUS_DISAGREE / STATUS_SINGLE, spelled out here
#: because family_registry imports this module and polish imports
#: family_registry (a test pins them equal).
MODELLED_STATUSES = frozenset({"polished_agree", "polished_disagree", "polished_single"})


@dataclass(frozen=True)
class SplitLocusParams:
    enabled: bool = True
    min_core_identity: float = 95.0
    max_edge_bp: int = 500
    min_flank_genes: int = 2
    flank_identity_gap: float = 10.0


#: Curator's ruling 2026-09-27; general (all families) unless a roster
#: `split_locus:` entry overrides it.
DEFAULT_SPLIT_LOCUS = SplitLocusParams()


def split_locus_params(locus: dict) -> SplitLocusParams:
    """Read a roster locus entry's optional `split_locus:` override.

    `split_locus: false` turns the rule off; a mapping overrides any subset of
    the parameters; absent means the defaults."""
    raw = locus.get("split_locus")
    if raw is None:
        return DEFAULT_SPLIT_LOCUS
    if raw is False:
        return SplitLocusParams(enabled=False)
    if raw is True:
        return DEFAULT_SPLIT_LOCUS
    known = {f.name for f in fields(SplitLocusParams)}
    return SplitLocusParams(**{k: v for k, v in raw.items() if k in known})


def _is_modelled(evidence) -> bool:
    return evidence.status in MODELLED_STATUSES or evidence.method == "diamond_proteome"


def _edge_distance(start: int, end: int, contig_length: int | None) -> int | None:
    if contig_length is None:
        return None
    return min(start - 1, contig_length - end)


def _overlaps(contig, start, end, spans) -> bool:
    return any(c == contig and start <= e and end >= s for c, s, e in spans)


def _strength(hit) -> tuple[float, float]:
    return (hit.bitscore if hit.bitscore is not None else -1.0, hit.identity)


def evaluate_split_locus(
    candidate,
    genome_hits,
    *,
    family_roles: dict[str, str],
    contig_lengths: dict[str, int],
    reported_spans: list[tuple[str, int, int]],
    params: SplitLocusParams = DEFAULT_SPLIT_LOCUS,
) -> dict | None:
    """The `split_locus` report block when `candidate` qualifies, else None.

    `candidate` is a built (withheld) DetectionResult for one cluster;
    `genome_hits` are this family's live hits across the whole genome;
    `family_roles` maps each roster gene to its role; `reported_spans` are
    (contig, start, end) of every locus already reported in this genome."""
    if not params.enabled:
        return None
    cores = [
        e for e in candidate.gene_evidence
        if e.role == "core_MAT" and _is_modelled(e)
        and e.identity is not None and e.identity >= params.min_core_identity
    ]
    if not cores:
        return None
    core = max(cores, key=lambda e: e.identity)
    edge = _edge_distance(core.start, core.end, contig_lengths.get(core.contig))
    if edge is None or edge > params.max_edge_bp:
        return None

    # Paralog guard: this must be the genome's strongest hit for the gene.
    same_gene = [h for h in genome_hits if h.gene_name == core.gene_name]
    if same_gene:
        best = max(same_gene, key=_strength)
        if not (best.contig == core.contig and best.start <= core.end and best.end >= core.start):
            return None

    flanks = []
    for gene, role in sorted(family_roles.items()):
        if not role.startswith("flanking"):
            continue
        hits = [h for h in genome_hits if h.gene_name == gene]
        if not hits:
            continue
        best = max(hits, key=_strength)
        if best.identity < core.identity - params.flank_identity_gap:
            continue
        if _overlaps(best.contig, best.start, best.end, reported_spans):
            continue
        flanks.append({"gene": gene, "contig": best.contig, "start": best.start,
                       "end": best.end, "identity": round(best.identity, 3)})
    if len(flanks) < params.min_flank_genes:
        return None
    if not any(f["contig"] != core.contig for f in flanks):
        return None
    return {
        "core_gene": core.gene_name,
        "core_contig": core.contig,
        "core_start": core.start,
        "core_end": core.end,
        "core_identity": round(core.identity, 3),
        "edge_distance": edge,
        "flanks": flanks,
        "contigs": sorted({core.contig} | {f["contig"] for f in flanks}),
    }
