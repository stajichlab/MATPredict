"""Report once a locus that several families called at the same place.

Curator's ruling 2026-09-27 (results/2026-09-27_caax_precursor/NOTE.md): the
Schizophyllum commune B locus came back three times over one span -- from the
generic Basidiomycota:PR family (receptor + CAAX-scan precursor) and from the
curated Balpha and Bbeta families. It is one physical locus, so it is reported
once, keeping every family's evidence.

What may merge is curation data, not code: a roster locus opts in with
`merge_group: <name>` and, if it is a broad catch-all for that locus type,
`merge_generic: true`. Only calls of the SAME group merge, so HD (group A)
never merges with PR (group B), and a family with no group is never touched.
Two calls merge only when they sit on the same contig, their spans overlap by
at least `MIN_SPAN_OVERLAP` of the shorter span, and their idiomorph labels do
not conflict (`undetermined` is compatible with anything). Genuine double
calls at separate loci -- e.g. the Cercospora kikuchii MAT1-1 and MAT1-2 calls
on two contigs -- therefore stay separate.

The merged call is the most specific member (non-generic first, then higher
confidence, then more genes): its family, and its idiomorph unless that is
`undetermined` and another member has a label. Its span, gene evidence,
genes found, reference records and segments are the union of the members';
its confidence is the best member confidence when all members carry the same
idiomorph, else the primary's; it is unverified if any member is (review
finding F7, 2026-09-28); `merged_from` records every member's family,
idiomorph, confidence and span.
"""
from __future__ import annotations

from dataclasses import replace

#: Two calls overlap when the shared span is at least this fraction of the
#: shorter one. The B-locus families report the identical cluster span, so any
#: value up to 1.0 merges them; 0.5 still keeps apart two loci that merely
#: abut or share a flank.
MIN_SPAN_OVERLAP = 0.5

_UNDETERMINED = {"undetermined", "", None}
_CONF_RANK = {"low": 0, "medium": 1, "high": 2}


def _key(family_key) -> str:
    return f"{family_key.phylum}:{family_key.locus_name}"


def _overlaps(a, b, min_overlap: float) -> bool:
    if a.contig != b.contig:
        return False
    shared = min(a.end, b.end) - max(a.start, b.start) + 1
    if shared <= 0:
        return False
    shorter = min(a.end - a.start, b.end - b.start) + 1
    return shared / shorter >= min_overlap


def _compatible(a, b) -> bool:
    return (a.idiomorph in _UNDETERMINED or b.idiomorph in _UNDETERMINED
            or a.idiomorph == b.idiomorph)


def _components(members: list[int], results: list, min_overlap: float) -> list[list[int]]:
    parent = {i: i for i in members}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for x, i in enumerate(members):
        for j in members[x + 1:]:
            if _overlaps(results[i], results[j], min_overlap) and _compatible(results[i], results[j]):
                parent[find(i)] = find(j)
    groups: dict[int, list[int]] = {}
    for i in members:
        groups.setdefault(find(i), []).append(i)
    return [sorted(g) for g in groups.values()]


def _merge(calls: list, groups: dict) -> object:
    def rank(r):
        generic = groups[r.family_key][1]
        return (generic, -_CONF_RANK.get(r.confidence, -1), -len(set(r.genes_found)), _key(r.family_key))

    primary = min(calls, key=rank)
    idiomorph = primary.idiomorph
    if idiomorph in _UNDETERMINED:
        labelled = [r for r in sorted(calls, key=rank) if r.idiomorph not in _UNDETERMINED]
        if labelled:
            idiomorph = labelled[0].idiomorph

    seen, evidence = set(), []
    for r in [primary] + [c for c in calls if c is not primary]:
        for e in r.gene_evidence:
            k = (e.gene_name, e.contig, e.start, e.end, e.method)
            if k not in seen:
                seen.add(k)
                evidence.append(e)

    def union(attr):
        out = []
        for r in [primary] + [c for c in calls if c is not primary]:
            for v in getattr(r, attr):
                if v not in out:
                    out.append(v)
        return out

    found = union("genes_found")
    # Review finding F7 (2026-09-28, results/2026-09-28_fable_review/): the
    # merged call is as cautious as its most cautious member. Verification:
    # unverified if ANY member is, reasons combined. Confidence: the best
    # member's only when every member carries the SAME idiomorph label (all
    # undetermined counts as the same); otherwise -- labels merely compatible,
    # one undetermined -- the primary's own confidence stands.
    if len({r.idiomorph for r in calls}) == 1:
        confidence = max((r.confidence for r in calls), key=lambda c: _CONF_RANK.get(c, -1))
    else:
        confidence = primary.confidence
    unverified = [r.verification for r in [primary] + [c for c in calls if c is not primary]
                  if r.verification and r.verification.get("status") == "unverified"]
    if not unverified:
        verification = primary.verification
    elif len(unverified) == 1:
        verification = unverified[0]
    else:
        reasons = []
        for v in unverified:
            if v.get("reason") and v["reason"] not in reasons:
                reasons.append(v["reason"])
        verification = {"status": "unverified", "reason": "; ".join(reasons),
                        "merged_from": unverified}
    return replace(
        primary,
        start=min(r.start for r in calls),
        end=max(r.end for r in calls),
        confidence=confidence,
        verification=verification,
        idiomorph=idiomorph,
        gene_evidence=evidence,
        genes_found=found,
        genes_missing=[g for g in primary.genes_missing if g not in found],
        reference_records=union("reference_records"),
        segments=union("segments"),
        polished_genes=max(r.polished_genes for r in calls),
        merged_from=[
            {"family": _key(r.family_key), "idiomorph": r.idiomorph,
             "confidence": r.confidence, "start": r.start, "end": r.end,
             "genes_found": list(r.genes_found)}
            for r in sorted(calls, key=rank)
        ],
    )


def merge_overlapping(results: list, groups: dict, min_overlap: float = MIN_SPAN_OVERLAP) -> list:
    """`results` with same-locus calls of one merge group collapsed.

    `groups` maps a FamilyKey to `(merge_group, merge_generic)`; families
    missing from it never merge. The first-reported position of each merged
    locus is kept, so unmerged calls keep their order.
    """
    by_bucket: dict[tuple, list[int]] = {}
    for i, r in enumerate(results):
        g = groups.get(r.family_key)
        if g and g[0]:
            by_bucket.setdefault((g[0], r.contig), []).append(i)

    replace_at: dict[int, object] = {}
    drop: set[int] = set()
    for members in by_bucket.values():
        if len(members) < 2:
            continue
        for comp in _components(members, results, min_overlap):
            calls = [results[i] for i in comp]
            if len(comp) < 2 or len({c.family_key for c in calls}) < 2:
                continue
            merged = _merge(calls, groups)
            first = comp[0]
            replace_at[first] = merged
            drop.update(i for i in comp if i != first)

    return [replace_at.get(i, r) for i, r in enumerate(results) if i not in drop]
