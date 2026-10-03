"""Simulate a per-family polish cap on an UNCAPPED run, for two rank rules.

A call is lost when no kept cluster of its family overlaps it. Validated on
f7b9773: genes-first N=6 reproduces the real cap6 run's 4 lost calls exactly.
It ignores knock-on effects between clusters, so treat it as an estimate.
usage: rank_simulation.py <dir of extracted per-genome reports> [...]
"""
import json, os, sys, collections, yaml
from yaml import CSafeLoader as L
def ov(a, b): return a["fam"] == b.get("family") and a["contig"] == b["contig"] and not (a["end"] < b["start"] or a["start"] > b["end"])
RANKS = {"genes-first": lambda c: (-c["genes"], -c["ident"], -c["hits"]),
         "identity-first": lambda c: (-c["ident"], -c["genes"], -c["hits"])}
for D in sys.argv[1:]:
    res = collections.defaultdict(collections.Counter); lost_by = collections.defaultdict(list); ng = 0
    for g in sorted(os.listdir(D)):
        rp = f"{D}/{g}/detection_report.yaml"
        if not os.path.exists(rp): continue
        ng += 1
        det = (yaml.load(open(rp), Loader=L) or {}).get("detected") or []
        ed = f"{D}/{g}/evidence_diagnostics.jsonl"
        adm = [dict(fam=e["family"], contig=e["contig"], start=e["cluster_start"], end=e["cluster_end"],
                    genes=e["gene_count"], ident=e["best_identity"], hits=e["hit_count"])
               for e in (map(json.loads, open(ed)) if os.path.exists(ed) else []) if e.get("kind") == "evidence" and e.get("admitted")]
        for rname, key in RANKS.items():
            for N in (4, 6, 8, 10):
                kept = []
                for f in {c["fam"] for c in adm}: kept += sorted([c for c in adm if c["fam"] == f], key=key)[:N]
                lost = [x for x in det if not any(ov(c, x) for c in kept)]
                r = res[(rname, N)]; r["lost"] += len(lost); r["calls"] += len(det)
                r["work"] += sum(c["genes"] for c in adm); r["kept"] += sum(c["genes"] for c in kept)
                lost_by[(rname, N)] += [(g[:15], x.get("family"), x["idiomorph"], x["confidence"], x.get("locus_class")) for x in lost]
    print(f"== {D}: {ng} genomes")
    for (rname, N), r in sorted(res.items()):
        print(f"  {rname:15s} N={N:2d}: calls lost {r['lost']:3d}/{r['calls']}  gene-load kept {100*r['kept']/max(r['work'],1):.0f}%")
    for k in (("genes-first", 6), ("identity-first", 6)):
        print("  ", k, lost_by[k][:30])
