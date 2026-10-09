#!/usr/bin/env python3
"""How often is a called locus held together by hits other than its own genes?

A cluster is chained from hits of every family with a maximum gap (family `max_cluster_gap_bp`, default 25 kb). If a called locus's
own gene evidence on its contig has an internal gap larger than that gap, something else (another family's hit, or an unreported
own-family hit) must have bridged it. Counted over saved campaign reports (called loci with >= 2 genes on the locus contig).
usage: bridging.py DB_ROOT OUT_TSV TARBALL [TARBALL ...]"""
import subprocess, sys, tarfile, io
from collections import Counter
import yaml
sys.path.insert(0, sys.argv[0].rsplit("/", 1)[0])
from MATPredict.detect.family_registry import load_all_families

db, out, tars = sys.argv[1], sys.argv[2], sys.argv[3:]
gap = {f"{f.key.phylum}:{f.key.locus_name}": f.max_cluster_gap_bp for f in load_all_families(__import__("pathlib").Path(db))}
rows = []
for tar in tars:
    name = tar.split("/")[-2]
    proc = subprocess.Popen(["zstd", "-dc", tar], stdout=subprocess.PIPE)
    with tarfile.open(fileobj=proc.stdout, mode="r|") as tf:
        for m in tf:
            if not m.name.endswith("detection_report.yaml"): continue
            y = yaml.safe_load(tf.extractfile(m).read()); g = m.name.split("/")[-2]
            for x in y.get("detected") or []:
                if x.get("fragmented"): continue
                genes = sorted((e["start"], e["end"]) for e in (x.get("gene_evidence") or []) if e.get("contig") == x["contig"])
                if len(genes) < 2: continue
                fam = x["family"]; lim = gap.get(fam, 25000)
                internal = max(b[0] - a[1] - 1 for a, b in zip(genes, genes[1:]))
                # genes can overlap or nest; use the running maximum end
                run, worst = genes[0][1], -10**9
                for s, e in genes[1:]:
                    worst = max(worst, s - run - 1); run = max(run, e)
                rows.append((name, g, fam, x["contig"], x["start"], x["end"], len(genes), worst, lim, worst > lim, x["confidence"], bool(x.get("merged_from"))))
with open(out, "w") as o:
    o.write("campaign\tgenome\tfamily\tcontig\tstart\tend\tn_genes\tmax_internal_gap\tfamily_gap\tbridged\tconfidence\tmerged\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
c = Counter((r[0], r[9]) for r in rows)
for camp in sorted({r[0] for r in rows}):
    n = c[(camp, True)] + c[(camp, False)]
    print(f"{camp}: called loci with >=2 genes {n}, internal gap > family gap: {c[(camp, True)]} ({100*c[(camp, True)]/max(n,1):.2f}%)")
print("total", len(rows), "bridged", sum(1 for r in rows if r[9]), "| bridged and merged_from set:", sum(1 for r in rows if r[9] and r[11]), "| bridged, not merged:", sum(1 for r in rows if r[9] and not r[11]))
for camp in sorted({r[0] for r in rows}):
    b = [r for r in rows if r[0] == camp and r[9] and not r[11]]
    print(f"{camp}: bridged and not merged {len(b)}; by family:", Counter(r[2] for r in b).most_common(5), "; median internal gap kb", sorted(r[7] for r in b)[len(b)//2]/1000 if b else None)
