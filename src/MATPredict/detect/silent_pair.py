"""`silent_pair`: split the merged S. pombe silent-cassette call into `mat2` (P) and `mat3` (M) (curator ruling 2026-11-05, #254).

*S. pombe* carries the expressed `mat1` locus (P or M) and two silent donor cassettes, `mat2`-P and `mat3`-M, about 11 kb apart. The
curated records are cassette-only (the mat1-P and mat2-P sequences are identical), so identity cannot tell the expressed locus from a silent
one; only the arrangement can. The detector reports the silent pair as ONE `mat1` call whose genes hold both a P cluster and an M cluster
(eight genomes of the v0.8.0 run). This rule splits such a call into two reported loci, `mat2` P and `mat3` M, when the two clusters lie
MIN_GAP_BP to MAX_GAP_BP apart on one contig, and marks each `cassette` with role `silent_pair`. A lone cassette elsewhere stays `mat1`
(the expressed locus). A call whose two clusters are farther apart than MAX_GAP_BP (GCA_006386235.1: 35 kb) is left as it is and marked
`arrangement_ambiguous`. The rule applies only where the species-scoped `mat2` and `mat3` families are in the genome's lineage, so no other
family or species is touched. Pure: the pipeline passes the calls, the families and the lineage.
"""
from __future__ import annotations

from dataclasses import replace

from MATPredict.detect.family_registry import FamilyKey

PHYLUM = "Ascomycota"
SOURCE, P_FAMILY, M_FAMILY = "mat1", "mat2", "mat3"
MIN_GAP_BP, MAX_GAP_BP = 5_000, 20_000
ROLE_PAIR, ROLE_AMBIGUOUS = "silent_pair", "arrangement_ambiguous"


def _scoped(family, lineage) -> bool:
    return family is not None and bool(set(family.taxonomic_scope) & set(lineage))


def _cluster(evidence):
    return (min(e.start for e in evidence), max(e.end for e in evidence)) if evidence else None


def split_silent_pairs(results: list, families: list, lineage: set[int]) -> list:
    """`results` with each merged mat1 call that holds a P and an M cluster 5-20 kb apart replaced by a `mat2` P and a `mat3` M call."""
    by_key = {f.key: f for f in families}
    src_key, p_key, m_key = (FamilyKey(PHYLUM, n) for n in (SOURCE, P_FAMILY, M_FAMILY))
    src, p_fam, m_fam = by_key.get(src_key), by_key.get(p_key), by_key.get(m_key)
    if not (_scoped(p_fam, lineage) and _scoped(m_fam, lineage) and src is not None):
        return list(results)
    side = {g["name"]: set(g.get("present_in_idiomorphs") or ()) for g in src.genes if g["role"] == "core_MAT"}
    out = []
    for r in results:
        if r.family_key != src_key:
            out.append(r)
            continue
        own = [e for e in r.gene_evidence if e.role == "core_MAT" and e.contig == r.contig and not getattr(e, "lost_to", None)]
        p_ev = [e for e in own if side.get(e.gene_name) == {"P"}]
        m_ev = [e for e in own if side.get(e.gene_name) == {"M"}]
        if not (p_ev and m_ev):
            out.append(r)
            continue
        (ps, pe), (ms, me) = _cluster(p_ev), _cluster(m_ev)
        gap = max(ps, ms) - min(pe, me)
        if not MIN_GAP_BP <= gap <= MAX_GAP_BP:
            out.append(replace(r, cassette={"idiomorph": r.idiomorph, "role": ROLE_AMBIGUOUS, "pair_distance_bp": gap,
                                            "note": f"a P and an M cluster {gap:,} bp apart (the silent pair is {MIN_GAP_BP:,}-{MAX_GAP_BP:,} bp); "
                                                    "not split, which cassette is expressed is not determined"}))
            continue
        for idio, key, fam, ev, other in (("P", p_key, p_fam, p_ev, "M"), ("M", m_key, m_fam, m_ev, "P")):
            lo, hi = _cluster(ev)
            renamed = [replace(e, gene_name=e.gene_name.replace("mat1-", "mat", 1)) for e in ev]
            found = sorted({e.gene_name for e in renamed})
            core = [g["name"] for g in fam.genes if g["role"] == "core_MAT" and not g.get("optional")]
            out.append(replace(
                r, family_key=key, start=lo, end=hi, idiomorph=idio, ambiguous_with=[], genes_found=found,
                genes_missing=[n for n in core if n not in found], gene_evidence=renamed,
                segments=[replace(s, start=lo, end=hi, contig_edge_distance=None) for s in r.segments[:1]] if r.segments else [],
                polished_genes=sum(1 for e in ev if str(e.status).startswith("polished")), supported_span=None,
                idiomorph_candidates=[c for c in (r.idiomorph_candidates or []) if c.get("idiomorph") == idio],
                idiomorph_class=(getattr(fam, "idiomorph_classes", {}) or {}).get(idio),
                cassette={"idiomorph": idio, "partner_idiomorphs": [other], "cross_family": True,
                          "partner_families": [f"{PHYLUM}:{m_key.locus_name if idio == 'P' else p_key.locus_name}"],
                          "role": ROLE_PAIR, "pair_distance_bp": gap,
                          "note": "the silent donor cassette of one pair: a P cluster (mat2) and an M cluster (mat3) about 11 kb apart; "
                                  "the expressed mat1 locus is a different cassette elsewhere"}))
    return out
