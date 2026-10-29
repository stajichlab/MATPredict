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
   genome is kept. If neither or both families are scoped to it, the call whose best core gene has the
   higher identity is kept (curator's decision, 2026-10-29: the confidence tier is not comparable across
   families, so it no longer decides). When the best identities are within `NEAR_TIE_POINTS` (3), the
   call from the family that declares no `confirm_outside_scope` rule is kept, else the better score
   (locus class, tier, polished genes, genes found, family name), and the kept call is flagged
   `family_unconfirmed` (reason `cross_call_near_tie`). The other call is not lost: it is recorded
   on the kept call as `cross_call`.
2. **Family status.** A call from a family that declares `confirm_outside_scope` (MTL) in a genome
   outside its scope is `family_unconfirmed` (reason `<family>_outside_scope`). It becomes
   `family_confirmed` when at least `min_flank_genes` of the declared flank genes (PAP1, OBP1, PIK1)
   were found within `window_bp` of a core gene of the call. The flag is separate from `confidence`,
   which is never changed. A genome with no lineage is never flagged: nothing says it is outside the scope.

3. **Cassettes.** A genome with two non-overlapping calls of a `cassette_report` family (MATog) of different idiomorphs
   has each call indicated as `cassette` (role `cassette_undetermined`): the MATa and MATalpha cassettes of a pseudohomothallic
   Ogataea. Which one is expressed is not determined, and the two calls are never merged.

Everything here is pure: the pipeline passes the calls, the families and the genome's lineage.
"""
from __future__ import annotations

from dataclasses import replace

POLICY_OFF = "off"
POLICY_ON = "on"
CROSS_FAMILY_POLICIES = (POLICY_OFF, POLICY_ON)
DEFAULT_CROSS_FAMILY_POLICY = POLICY_OFF

#: Best-core-identity difference (percentage points) at or below which two cross-called families are a near tie.
NEAR_TIE_POINTS = 3.0

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


def _best_core_identity(r) -> float:
    ids = [e.identity for e in r.gene_evidence if e.role == "core_MAT" and e.lost_to is None and e.identity is not None]
    return max(ids) if ids else 0.0


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

    def needs_confirmation(r) -> bool:
        family = by_key.get(r.family_key)
        return bool(getattr(family, "confirm_outside_scope", None))

    drop: set[int] = set()
    kept: dict[int, list[dict]] = {}
    for group in _groups(results):
        scoped = [i for i in group if in_scope(results[i])]
        near_tie = False
        if scoped:
            best = min(scoped, key=lambda i: _score(results[i]))
        else:
            top = max(_best_core_identity(results[i]) for i in group)
            close = [i for i in group if top - _best_core_identity(results[i]) <= NEAR_TIE_POINTS]
            near_tie = len(close) > 1
            pool = [i for i in close if not needs_confirmation(results[i])] or close
            best = min(pool, key=lambda i: _score(results[i])) if near_tie else close[0]
        notes = []
        for i in sorted(group):
            if i != best:
                note = _note(results[i])
                note["best_core_identity"] = round(_best_core_identity(results[i]), 2)
                if near_tie:
                    note["near_tie"] = True
                notes.append(note)
        kept[best] = notes
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
        if any(c.get("near_tie") for c in r.cross_call) and not (rule and lineage and not lineage & set(family.taxonomic_scope)):
            out.append(replace(r, family_status={"status": STATUS_UNCONFIRMED, "reason": "cross_call_near_tie"}))
            continue
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


CASSETTE_ROLE = "cassette_undetermined"


def flag_cassettes(results: list, families: list) -> list:
    """`results` with `cassette` set on the calls of a genome that carries two or more non-overlapping calls
    of one `cassette_report` family with different idiomorphs (the Ogataea MATa and MATalpha cassettes).

    Which cassette is the expressed locus is not determined, so every such call gets the role
    `cassette_undetermined`. A genome with a single call of the family, or two calls of one idiomorph, is not touched."""
    by_key = {f.key: f for f in families}
    out = list(results)
    for key in {r.family_key for r in results}:
        family = by_key.get(key)
        if not getattr(family, "cassette_report", False):
            continue
        idx = [i for i, r in enumerate(results) if r.family_key == key]
        idiomorphs = {results[i].idiomorph for i in idx}
        if len(idx) < 2 or len(idiomorphs) < 2:
            continue
        for i in idx:
            r = results[i]
            out[i] = replace(r, cassette={
                "idiomorph": r.idiomorph,
                "partner_idiomorphs": sorted(idiomorphs - {r.idiomorph}),
                "role": CASSETTE_ROLE,
                "note": "two cassettes of different idiomorphs in this genome; which one is expressed is not determined",
            })
    return out


def apply(results: list, families: list, lineage: set[int], policy: str) -> list:
    """The whole rule; `off` returns the calls untouched."""
    if check_policy(policy) == POLICY_OFF:
        return list(results)
    kept = resolve_cross_calls(results, families, lineage)
    return flag_cassettes(flag_family_status(kept, families, lineage), families)
