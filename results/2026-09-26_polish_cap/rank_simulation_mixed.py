"""Replay a per-family polish cap on UNCAPPED runs: genes-first, identity-first
and the MIXED rank the curator asked to test on 2026-09-26.

mixed = the top N-1 clusters by genes-first, plus the single highest-identity
cluster if it is not already among them (else the N-th genes-first cluster),
so it keeps exactly N clusters, like the other two ranks. mixed-4+2 (not
curator-ruled; measured as an alternative) = top N-2 by genes-first plus the
top 2 by identity, topped up to N by genes-first.

Same loss rule as rank_simulation.py (validated on f7b9773: genes-first N=6
reproduced the real cap6 run's 4 lost calls): a call is lost when no kept
cluster of its family overlaps it. Knock-on effects between clusters are
ignored, so treat the counts as estimates.

A call is also listed when one rank keeps it and another loses it, which is
the adoption test: mixed is adopted only if it loses nothing that either
single rank keeps.

usage: rank_simulation_mixed.py [--live] <runs dir> [...]
--live: rank on live distinct genes, as the real cap does (added after the
plain replay failed to reproduce the real Dothideomycetes losses).
"""
import collections, json, os, sys
import yaml
from yaml import CSafeLoader as L


def ov(a, b):
    return a["fam"] == b.get("family") and a["contig"] == b["contig"] and not (
        a["end"] < b["start"] or a["start"] > b["end"])


GF = lambda c: (-c["genes"], -c["ident"], -c["hits"])
IF = lambda c: (-c["ident"], -c["genes"], -c["hits"])


def keep(members, rank, n):
    if rank == "genes-first":
        return sorted(members, key=GF)[:n]
    if rank == "identity-first":
        return sorted(members, key=IF)[:n]
    k = 2 if rank == "mixed-4+2" else 1
    gf = sorted(members, key=GF)
    kept = gf[:n - k]
    for c in sorted(members, key=IF):
        if len(kept) >= n or sum(1 for x in kept if x in sorted(members, key=IF)[:k]) >= k:
            break
        if c not in kept:
            kept.append(c)
    for c in gf:                      # top up to n if the identity picks were already in
        if len(kept) >= n:
            break
        if c not in kept:
            kept.append(c)
    return kept


RANKS = ("genes-first", "identity-first", "mixed", "mixed-4+2")
NS = (6,)
LIVE = "--live" in sys.argv
for D in [a for a in sys.argv[1:] if not a.startswith("--")]:
    res = collections.defaultdict(collections.Counter)
    lost_by = collections.defaultdict(list)
    ng = 0
    for g in sorted(os.listdir(D)):
        rp = f"{D}/{g}/detection_report.yaml"
        if not os.path.exists(rp):
            continue
        ng += 1
        det = (yaml.load(open(rp), Loader=L) or {}).get("detected") or []
        ed = f"{D}/{g}/evidence_diagnostics.jsonl"
        rows = list(map(json.loads, open(ed))) if os.path.exists(ed) else []
        # The real cap ranks on LIVE distinct genes (superseded cross-hits
        # excluded, `_polish_rank`); the evidence row's gene_count includes
        # them. Pre-polish resolution events are the ones written before the
        # first evidence row; subtract their distinct losers per (family,
        # contig). Approximate when one contig holds several clusters.
        first_ev = next((i for i, e in enumerate(rows) if e.get("kind") == "evidence"), len(rows))
        losers = collections.defaultdict(set)
        for e in rows[:first_ev]:
            if e.get("kind") == "idiomorph_resolution":
                losers[(e["family"], e["contig"])].add(e["loser"])
        adm = [dict(fam=e["family"], contig=e["contig"], start=e["cluster_start"], end=e["cluster_end"],
                    genes=(e["gene_count"] - len(losers[(e["family"], e["contig"])]) if LIVE else e["gene_count"]),
                    ident=e["best_identity"], hits=e["hit_count"])
               for e in rows if e.get("kind") == "evidence" and e.get("admitted")]
        for N in NS:
            kept_by = {}
            for rname in RANKS:
                kept = []
                for f in {c["fam"] for c in adm}:
                    kept += keep([c for c in adm if c["fam"] == f], rname, N)
                kept_by[rname] = kept
                lost = [x for x in det if not any(ov(c, x) for c in kept)]
                r = res[(rname, N)]
                r["lost"] += len(lost)
                r["calls"] += len(det)
                r["work"] += sum(c["genes"] for c in adm)
                r["kept"] += sum(c["genes"] for c in kept)
                lost_by[(rname, N)] += [
                    (g, x.get("family"), x["contig"], x["start"], x["idiomorph"], x["confidence"])
                    for x in lost]
    print(f"== {D}: {ng} genomes")
    for (rname, N), r in sorted(res.items()):
        print(f"  {rname:15s} N={N:2d}: calls lost {r['lost']:3d}/{r['calls']}  "
              f"gene-load kept {100 * r['kept'] / max(r['work'], 1):.0f}%")
    for N in NS:
        sets = {k: set(lost_by[(k, N)]) for k in RANKS}
        kept_elsewhere = sets["mixed"] - (sets["genes-first"] & sets["identity-first"])
        print(f"  N={N}: mixed loses {len(sets['mixed'])}; of those, kept by a single rank: "
              f"{len(kept_elsewhere)}")
        for k in RANKS:
            for x in sorted(sets[k]):
                print(f"    lost[{k}] {x}")
