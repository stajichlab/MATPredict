#!/usr/bin/env python3
"""Ascomycota and Mucoromycota regression panels (110 genomes): does supported_span behave where there is no B-locus clustering?
Runs: base/ (current main, no new fields), floor33/ and floor39/ (supported_span code). usage: analyze_other.py (in this folder)"""
import glob, yaml
from collections import Counter

def load(d):
    return {f.split("/")[-2]: yaml.safe_load(open(f)) for f in glob.glob(f"{d}/runs/*/detection_report.yaml")}

def callkey(x):
    return (x["family"], x["contig"], x["start"], x["end"], x["confidence"], x["idiomorph"], x.get("locus_class"),
            tuple(sorted(x.get("genes_found") or [])), tuple(sorted(str(m) for m in (x.get("merged_from") or []))))

def sig(r):
    return ([callkey(x) for x in r.get("detected") or []],
            sorted((s["family"], s["contig"], s["start"], s["end"], s.get("withheld_reason"), tuple(s.get("genes_found") or [])) for s in r.get("suppressed_loci") or []))

import os
NAMES = [n for n in ("base", "floor33", "floor39") if glob.glob(f"{n}/runs/*/detection_report.yaml")]
D = {n: load(n) for n in NAMES}
print("reports:", {n: len(v) for n, v in D.items()})
phylum = lambda fam: fam.split(":")[0]
print("\n1. CALL IDENTITY against current main (calls, cluster spans, genes, merged_from, withheld loci)")
for n in [x for x in NAMES if x != "base"]:
    common = [g for g in D[n] if g in D["base"]]
    bad = [g for g in common if sig(D[n][g]) != sig(D["base"][g])]
    print(f"   {n}: {len(common)} common genomes, {len(bad)} with any difference {bad[:4]}")

for n in [x for x in NAMES if x != "base"]:
    print(f"\n--- {n} ---")
    for ph in ("Ascomycota", "Mucoromycota"):
        called = [(g, x) for g, r in D[n].items() for x in r.get("detected") or [] if phylum(x["family"]) == ph]
        withheld = [(g, x) for g, r in D[n].items() for x in r.get("suppressed_loci") or [] if phylum(x["family"]) == ph]
        cb = [x["core_span"]["beyond_core_bp"] for _, x in called if x.get("core_span")]
        sb = [x["supported_span"]["beyond_supported_bp"] for _, x in called if x.get("supported_span")]
        inv = sum(1 for _, x in called if x.get("core_span") and x.get("supported_span") and not (x["supported_span"]["start"] <= x["core_span"]["start"] and x["core_span"]["end"] <= x["supported_span"]["end"]))
        out = sum(1 for _, x in called for e in x.get("gene_evidence") or [] if e.get("contig") == x["contig"] and x.get("supported_span") and not (x["supported_span"]["start"] <= e["start"] and e["end"] <= x["supported_span"]["end"]))
        ng = sum(1 for _, x in called for e in x.get("gene_evidence") or [] if e.get("contig") == x["contig"])
        print(f"{ph}: called {len(called)}; with core/supported span {len(cb)}/{len(sb)}; genes checked {ng}, outside supported span {out}; core not inside supported {inv}")
        for lab, v in (("beyond_core", cb), ("beyond_supported", sb)):
            if v:
                print(f"   {lab:17s} >0: {sum(1 for a in v if a>0):4d}  >=1kb: {sum(1 for a in v if a>=1000):4d}  >=10kb: {sum(1 for a in v if a>=10000):4d}  total {sum(v)/1000:8.1f} kb  max {max(v)/1000:.1f} kb")
        wc = [x["core_span"]["beyond_core_bp"] for _, x in withheld if x.get("core_span")]
        ws = [x["supported_span"]["beyond_supported_bp"] for _, x in withheld if x.get("supported_span")]
        print(f"   withheld loci {len(withheld)}: with core_span {len(wc)}, beyond_core >=10kb {sum(1 for a in wc if a>=10000)}; with supported_span {len(ws)}, beyond_supported >=10kb {sum(1 for a in ws if a>=10000)}")
    if n == "floor33":
        big = sorted(((x["supported_span"]["beyond_supported_bp"], g, x) for g, r in D[n].items() for x in r.get("detected") or []
                      if x.get("supported_span") and x["supported_span"]["beyond_supported_bp"] >= 10000), key=lambda t: -t[0])
        print("   called loci with >= 10 kb beyond the supported span (all phyla):", len(big))
        for b, g, x in big[:12]:
            print(f"      {g[:30]:30s} {x['family']:22s} {x['contig']} cluster {x['start']}-{x['end']} supported {x['supported_span']['start']}-{x['supported_span']['end']} beyond {b} genes {x['genes_found']}")
