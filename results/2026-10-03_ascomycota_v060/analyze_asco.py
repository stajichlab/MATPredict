#!/usr/bin/env python3
"""Aggregate the full Ascomycota v0.6.0 run.

Reads wave_*/runs/*/detection_report.yaml and wall_seconds; joins BFD samples.csv.
Writes genomes.tsv (one row per listed genome), loci.tsv (one row per call) and
by_class.tsv (genomes, reports, called, routing modes, wall median/p90/max).
Run with /usr/bin/python3.12.
"""
import csv, glob, os, statistics
from collections import Counter, defaultdict
import yaml
try:
    from yaml import CSafeLoader as L
except ImportError:
    from yaml import SafeLoader as L

H = os.path.dirname(os.path.abspath(__file__))
samp = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
listed = [l.split("\t")[0] for f in sorted(glob.glob(f"{H}/lists/wave_*.tsv")) for l in open(f) if l.strip()]
dirs = {os.path.basename(d): d for d in glob.glob(f"{H}/wave_*/runs/*")}
G, LO = [], []
for a in listed:
    s = samp.get(a, {}); d = dirs.get(a); rep = f"{d}/detection_report.yaml" if d else None
    row = dict(genome=a, species=s.get("SPECIES", ""), class_=s.get("CLASS", "") or "NA", order=s.get("ORDER", ""),
               family=s.get("FAMILY", ""), status="no_report", routing="", n_loci=0, families_called="", wall_s="")
    if d and os.path.exists(f"{d}/wall_seconds"):
        row["wall_s"] = open(f"{d}/wall_seconds").read().strip()
    if rep and os.path.exists(rep) and os.path.getsize(rep):
        doc = yaml.load(open(rep), Loader=L) or {}
        det = doc.get("detected") or []
        r = doc.get("routing")
        row.update(status="called" if det else "uncalled", n_loci=len(det),
                   routing=str(r.get("mode")) if isinstance(r, dict) else str(doc.get("routing_mode", "")),
                   families_called="|".join(sorted({x.get("family") for x in det})))
        for x in det:
            v = x.get("verification")
            LO.append(dict(genome=a, class_=row["class_"], order=row["order"], family_called=x.get("family"),
                           contig=x.get("contig"), start=x.get("start"), end=x.get("end"), idiomorph=x.get("idiomorph"),
                           confidence=x.get("confidence"), locus_class=x.get("locus_class"),
                           detection_pass=x.get("detection_pass"), genes_found="|".join(x.get("genes_found") or []),
                           verification=(v.get("status") if isinstance(v, dict) else (v or ""))))
    elif d and os.path.exists(f"{d}/stderr.log"):
        row["status"] = "failed"
    G.append(row)
for name, rows in (("genomes.tsv", G), ("loci.tsv", LO)):
    with open(f"{H}/{name}", "w") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
by = defaultdict(list)
for r in G: by[r["class_"]].append(r)
def q(v, p): v = sorted(v); return v[min(len(v) - 1, int(round((len(v) - 1) * p)))]
with open(f"{H}/by_class.tsv", "w") as fo:
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(["class", "genomes", "reports", "called", "pct_called", "routing", "wall_median_s", "wall_p90_s", "wall_max_s"])
    for c, rs in sorted(by.items(), key=lambda x: -len(x[1])):
        rep = [r for r in rs if r["status"] in ("called", "uncalled")]
        wall = [int(r["wall_s"]) for r in rep if r["wall_s"]]
        called = sum(r["status"] == "called" for r in rs)
        w.writerow([c, len(rs), len(rep), called, f"{100 * called / max(1, len(rep)):.1f}",
                    "; ".join(f"{k}:{n}" for k, n in Counter(r["routing"] for r in rep).most_common()),
                    int(statistics.median(wall)) if wall else "", q(wall, .9) if wall else "", max(wall) if wall else ""])
print(Counter(r["status"] for r in G), len(LO), "loci")
print("families:", Counter(l["family_called"] for l in LO).most_common(25))
print("verification:", Counter(l["verification"] for l in LO))
wall = [int(r["wall_s"]) for r in G if r["wall_s"] and r["status"] != "failed"]
print(f"wall median {statistics.median(wall):.0f} s, total {sum(wall)/3600:.0f} CPU-h")
