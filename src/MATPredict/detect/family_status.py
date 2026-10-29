"""Cross-family calls and the `family_unconfirmed` flag, resolved after searching.

Curator's decisions, 2026-10-29 (issues #188 and #189; evidence in
results/2026-10-29_yeast_routing/). Behind `--cross-family-policy {off,on}`
(default `off`, so a run without the flag is unchanged).

A genome routed by phylum fallback is searched with every Ascomycota family, so one yeast locus can
be called by two of them (v0.7.0: 181 genomes carried a MATsc and an MTL locus at the same place) and
an MTL call can land on a genome far outside the Serinales (22 MTL-only calls, 5 with the C. albicans
PAP1, OBP1, PIK1 neighbours, 17 without). With the policy `on`, after every call decision:

1. **Cross-calls.** Two families calling overlapping loci (overlap of 0.9 or more of the shorter
   locus, on one contig) are one locus. The call from the family whose `taxonomic_scope` contains the
   genome is kept; if neither or both families are scoped to it, the better-scoring call is kept
   (locus class, then confidence tier, then polished genes, then genes found, then the family name).
   The other call is not lost: it is recorded on the kept call as `cross_call`.
2. **Family status.** A call from a family that declares `confirm_outside_scope` (MTL) in a genome
   outside its scope is `family_unconfirmed` (reason `<family>_outside_scope`). It becomes
   `family_confirmed` when at least `min_flank_genes` of the declared flank genes (PAP1, OBP1, PIK1)
   were found within `window_bp` of a core gene of the call. The flag is separate from `confidence`,
   which is never changed. A genome with no lineage is never flagged: nothing says it is outside the scope.

Everything here is pure: the pipeline passes the calls, the families and the genome's lineage.
"""
from __future__ import annotations

from dataclasses import replace

POLICY_OFF = "off"
POLICY_ON = "on"
CROSS_FAMILY_POLICIES = (POLICY_OFF, POLICY_ON)
DEFAULT_CROSS_FAMILY_POLICY = POLICY_OFF

#: Fraction of the shorter locus that two calls must share to count as one locus.
MIN_OVERLAP = 0.9

STATUS_UNCONFIRMED = "family_unconfirmed"
STATUS_CONFIRMED = "family_confirmed"

_CLASS_RANK = {"mat_locus": 0, "partial_locus": 1, "idiomorph_gene_only": 2}
_TIER_RANK = {"high": 0, "medium": 1, "low": 2}


def check_policy(policy: str) -> str:
    if policy not in CROSS_FAMILY_POLICIES:
        raise ValueError(f"unknown cross-family policy {policy!r}; expected one of {CROSS_FAMILY_POLICIES}")
    return policy


def _overlap_fraction(a, b) -> float:
    """Overlap of two calls on one contig as a fraction of the shorter one (0 on different contigs)."""
    if a.contig != b.contig:
        return 0.0
    shared = min(a.end, b.end) - max(a.start, b.start) + 1
    shorter = min(a.end - a.start, b.end - b.start) + 1
    return max(0, shared) / shorter if shorter > 0 else 0.0


def _polished(r) -> int:
    return sum(1 for e in r.gene_evidence if str(e.status).startswith("polished") and e.lost_to is None)


def _score(r) -> tuple:
    """Smaller is better; the same order as `scripts/dedupe_family_calls.py`."""
    return (_CLASS_RANK.get(r.locus_class, 9), _TIER_RANK.get(r.confidence, 9), -_polished(r),
            -len(r.genes_found), r.family_key.phylum, r.family_key.locus_name)


def _groups(results: list) -> list[list[int]]:
    """Indices of calls that overlap one another across different families (union-find)."""
    parent = list(range(len(results)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(results)):
        for j in range(i + 1, len(results)):
            a, b = results[i], results[j]
            if a.family_key != b.family_key and _overlap_fraction(a, b) >= MIN_OVERLAP:
                parent[find(i)] = find(j)
    out: dict[int, list[int]] = {}
    for i in range(len(results)):
        out.setdefault(find(i), []).append(i)
    return [g for g in out.values() if len({results[i].family_key for i in g}) > 1]


def _note(r) -> dict:
    return {"family": f"{r.family_key.phylum}:{r.family_key.locus_name}", "idiomorph": r.idiomorph,
            "confidence": r.confidence, "locus_class": r.locus_class, "contig": r.contig,
            "start": r.start, "end": r.end, "genes_found": list(r.genes_found)}


def resolve_cross_calls(results: list, families: list, lineage: set[int]) -> list:
    """`results` with each cross-family group reduced to one call; the others ride on it as `cross_call`."""
    by_key = {f.key: f for f in families}

    def in_scope(r) -> bool:
        family = by_key.get(r.family_key)
        return bool(lineage and family and lineage & set(family.taxonomic_scope))

    drop: set[int] = set()
    kept: dict[int, list[dict]] = {}
    for group in _groups(results):
        best = min(group, key=lambda i: (0 if in_scope(results[i]) else 1, _score(results[i])))
        kept[best] = [_note(results[i]) for i in sorted(group) if i != best]
        drop.update(i for i in group if i != best)
    return [replace(r, cross_call=kept[i]) if i in kept else r for i, r in enumerate(results) if i not in drop]


def _gap(a, b) -> int:
    """Base pairs between two gene spans on one contig (0 when they touch or overlap)."""
    return max(0, max(a.start, b.start) - min(a.end, b.end) - 1)


def flank_support(r, flank_genes, window_bp: int) -> list[str]:
    """The declared flank genes that call `r` carries within `window_bp` of one of its core genes."""
    wanted = set(flank_genes)
    core = [e for e in r.gene_evidence if e.role == "core_MAT" and e.lost_to is None]
    found = {e.gene_name for e in r.gene_evidence
             if e.gene_name in wanted and e.lost_to is None
             and any(c.contig == e.contig and _gap(c, e) <= window_bp for c in core)}
    return sorted(found)


def flag_family_status(results: list, families: list, lineage: set[int]) -> list:
    """`results` with `family_status` set on every call whose family declares `confirm_outside_scope`
    and whose genome lies outside that family's scope."""
    by_key = {f.key: f for f in families}
    out = []
    for r in results:
        family = by_key.get(r.family_key)
        rule = getattr(family, "confirm_outside_scope", None)
        if not rule or not lineage or lineage & set(family.taxonomic_scope):
            out.append(r)
            continue
        genes = list(rule["flank_genes"])
        window = int(rule.get("window_bp", 15000))
        need = int(rule.get("min_flank_genes", 2))
        found = flank_support(r, genes, window)
        out.append(replace(r, family_status={
            "status": STATUS_CONFIRMED if len(found) >= need else STATUS_UNCONFIRMED,
            "reason": f"{family.key.locus_name.lower()}_outside_scope",
            "flank_genes_found": found,
            "flank_genes_declared": genes,
            "min_flank_genes": need,
            "window_bp": window,
        }))
    return out


def apply(results: list, families: list, lineage: set[int], policy: str) -> list:
    """The whole rule; `off` returns the calls untouched."""
    if check_policy(policy) == POLICY_OFF:
        return list(results)
    return flag_family_status(resolve_cross_calls(results, families, lineage), families, lineage)
