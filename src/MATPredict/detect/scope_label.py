"""Label a Basidiomycota MAT call outside its scope that rests on weak genes only (`--scope-label-policy`, #170).

Curator's decision, 2026-11-07 (issue #170): label as unverified, keep the call. Behind `--scope-label-policy {off,on}` (default `off`).

The `Basidiomycota:MAT` family (the Tremellomycetes-type locus) is called in other classes: in the v0.8.1 run 13 loci of Agaricomycetes,
Ustilaginomycetes, Microbotryomycetes and Exobasidiomycetes, 10 of them on STE3 plus PAN6 or RPL39 alone (*Ceratobasidium* 7, *Thecaphora* 2,
*Neolentinus* 1; 3 are `high`, one on PAN6 at 10.9% coverage and STE3 at 37.6%). PAN6 and RPL39 are conserved genes that sit near the
Tremellomycetes MAT locus but are not mating-type genes, and STE3 is a pheromone-receptor family member present in every basidiomycete.

With the policy `on`, such a call is marked with the existing `family_status` block: status `family_unconfirmed`, reason
`mat_outside_scope_weak_genes_only`. The call, its idiomorph, coordinates and tier are never changed (the same report-only pattern as the
MTL flag of `--cross-family-policy`). A genome with an unknown lineage is never flagged: nothing says it is outside the scope.
"""
from __future__ import annotations

from dataclasses import replace

POLICY_OFF = "off"
POLICY_ON = "on"
OUT_OF_SCOPE_LABEL_POLICIES = (POLICY_OFF, POLICY_ON)
DEFAULT_OUT_OF_SCOPE_LABEL_POLICY = POLICY_OFF

SCOPE_FAMILY = ("Basidiomycota", "MAT")
#: A call whose genes are all among these, in a genome outside the family's scope, is unverified.
WEAK_GENES = frozenset({"STE3", "PAN6", "RPL39"})
REASON = "mat_outside_scope_weak_genes_only"
STATUS_UNCONFIRMED = "family_unconfirmed"


def check_policy(policy: str) -> str:
    if policy not in OUT_OF_SCOPE_LABEL_POLICIES:
        raise ValueError(f"unknown out-of-scope-label policy {policy!r}; expected one of {OUT_OF_SCOPE_LABEL_POLICIES}")
    return policy


def apply(results: list, families: list, lineage: set[int], policy: str) -> list:
    """`results` with `family_status` set on each such call; `off` or an unknown lineage returns them untouched."""
    if check_policy(policy) == POLICY_OFF or not lineage:
        return list(results)
    by_key = {f.key: f for f in families}
    out = []
    for r in results:
        family = by_key.get((r.family_key))
        if (family is None or (r.family_key.phylum, r.family_key.locus_name) != SCOPE_FAMILY or getattr(r, "family_status", None)
                or not family.taxonomic_scope or lineage & set(family.taxonomic_scope)):
            out.append(r)
            continue
        genes = set(r.genes_found)
        if not genes or not genes <= WEAK_GENES:
            out.append(r)
            continue
        out.append(replace(r, family_status={
            "status": STATUS_UNCONFIRMED, "reason": REASON, "genes_found": sorted(genes), "weak_genes": sorted(WEAK_GENES),
            "family_taxonomic_scope": list(family.taxonomic_scope),
        }))
    return out
