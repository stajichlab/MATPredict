#!/usr/bin/env python
"""For loci the v0.6.0 replay predicted to be capped, show what the current-code `partial_weak` report says about the same locus."""
import sys, yaml
from pathlib import Path
L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
runs = Path(sys.argv[1])
for line in open(sys.argv[2]):
    g, fam, contig, start, bg, cov, ident, edge, nfl = line.rstrip("\n").split("\t")
    doc = yaml.load((runs / "partial_weak" / g / "detection_report.yaml").open(), Loader=L)
    m = [c for c in doc.get("detected") or [] if c["family"] == fam and c["contig"] == contig and abs(c["start"] - int(start)) < 10000]
    if not m:
        print(g[:28], fam[9:], "v060:", bg, cov, ident, "edge", edge, "-> NO MATCHING LOCUS now"); continue
    c = m[0]
    ge = [(e["gene"], e.get("role", "")[:4], e["identity"], e.get("coverage"), e["status"][-6:]) for e in c["gene_evidence"] if e["gene"] == bg]
    print(g[:28], fam[9:], "v060:", bg, cov, ident, "edge", edge, "fl", nfl, "-> now", c["confidence"], c["locus_class"], "cap" if c.get("confidence_cap") else "", ge[:1], "segedge", [s.get("contig_edge_distance") for s in c["segments"]])
