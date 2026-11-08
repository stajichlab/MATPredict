"""A `high` Basidiomycota call must rest on two distinct core regions or a polished flank gene (`--single-region-cap-policy`, #167).

Curator's decision, 2026-11-07 (issue #167): opt-in cap for every Basidiomycota family. Behind `--single-region-cap-policy {off,on}` (default `off`).

In the Basidiomycota HD locus the HD1 and HD2 references both hit the same homeodomain region (their models overlap by half or more in
1,063 of 2,529 HD loci), so one region can satisfy two core genes; bE fragments covering 10-13% of the protein have reached `high` the same
way. Measured on the v0.8.1 Basidiomycota run: 202 `high` loci rest on fewer than two distinct core regions with no polished flank gene (HD 196,
194 of them Agaricomycetes; bLocus 4; Aalpha 1; redHD 1), of 3,622 `high` Basidiomycota loci.

With the policy `on`, a `high` Basidiomycota call whose live core genes fall in fewer than two distinct regions, and which carries no polished flank
gene, is lowered to `medium`. Models on one contig that overlap by `MIN_OVERLAP` (half) of the shorter one are one region. Only the tier changes;
the idiomorph, the locus and the coordinates never do. The report's `single_region_cap` says what was found.
"""
from __future__ import annotations

from dataclasses import dataclass

POLICY_OFF = "off"
POLICY_ON = "on"
SINGLE_REGION_CAP_POLICIES = (POLICY_OFF, POLICY_ON)
DEFAULT_SINGLE_REGION_CAP_POLICY = POLICY_OFF

PHYLUM = "Basidiomycota"
CORE_ROLE = "core_MAT"
FLANK_ROLE_PREFIX = "flanking"
#: Two models on one contig overlapping by at least this fraction of the shorter are one region.
MIN_OVERLAP = 0.5


def check_policy(policy: str) -> str:
    if policy not in SINGLE_REGION_CAP_POLICIES:
        raise ValueError(f"unknown single-region-cap policy {policy!r}; expected one of {SINGLE_REGION_CAP_POLICIES}")
    return policy


@dataclass(frozen=True)
class GeneStat:
    """What the rule needs to know of one live gene (not `lost_to`)."""

    gene: str
    role: str | None
    status: str | None
    contig: str
    start: int
    end: int


def count_regions(core: list[GeneStat]) -> int:
    """Distinct regions among `core`: models on one contig that overlap by >= `MIN_OVERLAP` of the shorter are one."""
    regions: list[list] = []  # [contig, start, end]
    for g in core:
        for r in regions:
            if r[0] != g.contig:
                continue
            shared = min(r[2], g.end) - max(r[1], g.start) + 1
            shorter = min(r[2] - r[1], g.end - g.start) + 1
            if shared > 0 and shared / shorter >= MIN_OVERLAP:
                r[1], r[2] = min(r[1], g.start), max(r[2], g.end)
                break
        else:
            regions.append([g.contig, g.start, g.end])
    return len(regions)


def evaluate(confidence: str, phylum: str, stats: list[GeneStat], *, policy: str = POLICY_ON) -> dict | None:
    """None when the call keeps its tier; otherwise the note saying why it drops from `high` to `medium`."""
    if policy == POLICY_OFF or confidence != "high" or phylum != PHYLUM:
        return None
    core = [s for s in stats if s.role == CORE_ROLE]
    if not core:
        return None
    if any(str(s.role or "").startswith(FLANK_ROLE_PREFIX) and str(s.status or "").startswith("polished") for s in stats):
        return None
    regions = count_regions(core)
    if regions >= 2:
        return None
    return {"rule": "single_region_cap", "from": "high", "to": "medium", "core_genes": sorted({s.gene for s in core}),
            "distinct_core_regions": regions, "min_overlap_fraction": MIN_OVERLAP, "polished_flank_gene": False}
