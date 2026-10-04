#!/usr/bin/env python3
"""R. toruloides / Sporidiobolales HD: generic-HD call (2026-09-26, run-ad1f865)
versus redHD call (v0.6.0, run-7c7ed99), per genome.

For each Sporidiobolales genome in both runs: the old Basidiomycota:HD call
(contig, span, HD1/HD2/MIP1 identity, model status, HD1 model length) and the
new Basidiomycota:redHD call (contig, span, HD1/HD2 identity and best reference
record), plus the redPR contig. Writes per_genome.tsv; prints summaries.
"""
import csv, glob, os, statistics, collections
import yaml
try:
    from yaml import CSafeLoader as L
except ImportError:
    from yaml import SafeLoader as L
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
H = f"{R}/2026-10-04_rtoruloides_hd"
G = {r["genome"]: r for r in csv.DictReader(open(f"{R}/2026-10-03_basidiomycota_v060/genomes.tsv"), delimiter="\t")
     if r["order"] == "Sporidiobolales"}
def rep(run, g):
    f = glob.glob(f"{R}/{run}/*/runs/{g}/detection_report.yaml")
    return yaml.load(open(f[0]), Loader=L) if f else None
def ev(call, gene):
    es = [e for e in call.get("gene_evidence") or [] if e.get("gene") == gene]
    return max(es, key=lambda e: e.get("bitscore") or 0) if es else None
rows = []
for g, gr in sorted(G.items()):
    o, n = rep("2026-09-26_basidiomycota_full", g), rep("2026-10-03_basidiomycota_v060", g)
    row = dict(genome=g, species=gr["species"])
    if o:
        hd = [x for x in o.get("detected") or [] if x["family"] == "Basidiomycota:HD"]
        row["old_routing"] = o.get("routing_mode")
        if hd:
            x = hd[0]
            row.update(old_contig=x["contig"], old_start=x["start"], old_end=x["end"], old_conf=x.get("confidence"),
                       old_genes="|".join(x.get("genes_found") or []))
            for gene in ("HD1", "HD2", "MIP1"):
                e = ev(x, gene)
                if e:
                    row[f"old_{gene}_id"] = e.get("identity")
                    row[f"old_{gene}_status"] = e.get("status")
                    row[f"old_{gene}_span"] = e["end"] - e["start"] + 1
    if n:
        row["new_routing"] = n.get("routing_mode")
        for fam, key in (("Basidiomycota:redHD", "new"), ("Basidiomycota:redPR", "pr")):
            xs = [x for x in n.get("detected") or [] if x["family"] == fam]
            if xs:
                x = xs[0]
                row.update({f"{key}_contig": x["contig"], f"{key}_start": x["start"], f"{key}_end": x["end"],
                            f"{key}_conf": x.get("confidence"), f"{key}_idiomorph": x.get("idiomorph")})
                if key == "new":
                    for gene in ("HD1", "HD2"):
                        e = ev(x, gene)
                        if e:
                            row[f"new_{gene}_id"] = e.get("identity")
                            row[f"new_{gene}_ref"] = e.get("reference_record")
                            row[f"new_{gene}_status"] = e.get("status")
    if row.get("old_contig") and row.get("new_contig"):
        same = row["old_contig"] == row["new_contig"] and int(row["old_start"]) <= int(row["new_end"]) and int(row["new_start"]) <= int(row["old_end"])
        row["old_vs_new"] = "same_locus" if same else ("same_contig" if row["old_contig"] == row["new_contig"] else "different_contig")
    if row.get("old_contig") and row.get("pr_contig"):
        row["old_on_pr_contig"] = row["old_contig"] == row["pr_contig"]
    if row.get("new_contig") and row.get("pr_contig"):
        row["new_on_pr_contig"] = row["new_contig"] == row["pr_contig"]
    rows.append(row)
keys = []
for r in rows:
    for k in r:
        if k not in keys: keys.append(k)
with open(f"{H}/per_genome.tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=keys, delimiter="\t", lineterminator="\n", restval=""); w.writeheader(); w.writerows(rows)
C = collections.Counter
print("genomes", len(rows))
print("old HD called", sum(bool(r.get("old_contig")) for r in rows), "| new redHD called", sum(bool(r.get("new_contig")) for r in rows))
print("old vs new:", C(r.get("old_vs_new") for r in rows if r.get("old_contig")))
print("old HD on the redPR contig:", C(r.get("old_on_pr_contig") for r in rows if r.get("old_contig")))
print("new redHD on the redPR contig:", C(r.get("new_on_pr_contig") for r in rows if r.get("new_contig")))
def s(v):
    v = [float(x) for x in v if x not in (None, "")]
    return f"n={len(v)} median={statistics.median(v):.1f} min={min(v):.1f} max={max(v):.1f}" if v else "n=0"
for tag, sel in (("old (generic HD)", lambda r: r.get("old_contig")), ("new (redHD)", lambda r: r.get("new_contig"))):
    rr = [r for r in rows if sel(r)]
    pre = "old" if tag.startswith("old") else "new"
    print(tag, "HD1 identity", s(r.get(f"{pre}_HD1_id") for r in rr), "| HD2 identity", s(r.get(f"{pre}_HD2_id") for r in rr))
print("old HD1 model span bp", s(r.get("old_HD1_span") for r in rows))
print("old HD1 status", C(r.get("old_HD1_status") for r in rows if r.get("old_contig")), "HD2 status", C(r.get("old_HD2_status") for r in rows if r.get("old_contig")))
print("new HD1 status", C(r.get("new_HD1_status") for r in rows if r.get("new_contig")))
print("old routing", C(r.get("old_routing") for r in rows), "new routing", C(r.get("new_routing") for r in rows))
print("by species (n, old called, new called, same locus):")
sp = collections.defaultdict(lambda: [0, 0, 0, 0])
for r in rows:
    k = sp[r["species"]]; k[0] += 1; k[1] += bool(r.get("old_contig")); k[2] += bool(r.get("new_contig")); k[3] += r.get("old_vs_new") == "same_locus"
for k, v in sorted(sp.items(), key=lambda x: -x[1][0])[:12]:
    print("  ", k, v)
