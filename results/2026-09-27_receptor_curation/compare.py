"""Before/after comparison for the Polyporales/Russulales PR records.

before = run_f26c397 (no Polyporales/Russulales PR records, PR scope Agaricales)
after  = run_29818af (three tier-2 PR records, PR scope + Polyporales, Russulales)
Reports PR and HD calls per genome, whether each PR call overlaps a strict-CAAX
flagged STE3 locus from the positional scan, whether PR and HD share a contig,
receptor/precursor copies per call, and wall time.
"""
import csv
import glob
import os
import statistics
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
POS = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_pheromone_positional/uncurated_loci.tsv"
samples = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
flag = {}
for r in csv.DictReader(open(POS), delimiter="\t"):
    if r["flag_strict"] == "1":
        flag.setdefault(r["asm"], []).append((r["contig"], int(r["start"]), int(r["end"])))
panel = [l.split("\t")[0] for l in open(f"{HERE}/lists/panel.tsv")]


def load(arm):
    out = {}
    for a in panel:
        p = f"{HERE}/run_{arm}/runs/{a}/detection_report.yaml"
        if not os.path.exists(p):
            out[a] = None
            continue
        r = yaml.safe_load(open(p)) or {}
        w = f"{HERE}/run_{arm}/runs/{a}/wall_seconds"
        r["_wall"] = float(open(w).read().strip()) if os.path.exists(w) else None
        out[a] = r
    return out


def calls(r, fam):
    return [x for x in (r.get("detected") or []) if x["family"] == fam]


B, A = load("f26c397"), load("29818af")
rows = []
for a in panel:
    s = samples.get(a, {})
    rb, ra = B[a], A[a]
    if ra is None:
        rows.append(dict(asm=a, order=s.get("ORDER"), species=s.get("SPECIES"), status="no report"))
        continue
    pr_a = calls(ra, "Basidiomycota:PR")
    hd_a = calls(ra, "Basidiomycota:HD")
    pr_b = calls(rb, "Basidiomycota:PR") if rb else []
    hd_b = calls(rb, "Basidiomycota:HD") if rb else []
    fl = flag.get(a, [])
    at_flag = sum(any(c == x["contig"] and s0 <= x["end"] and e0 >= x["start"] for c, s0, e0 in fl) for x in pr_a)
    hd_contigs = {x["contig"] for x in hd_a}
    shared = sum(x["contig"] in hd_contigs for x in pr_a)
    ev = [e for x in pr_a for e in x.get("gene_evidence", [])]
    nrec = sum(e["gene"] == "pheromone_receptor" for e in ev)
    npre = sum(e["gene"] == "fungal_mating_type_pheromone" for e in ev)
    rows.append(dict(
        asm=a, order=s.get("ORDER"), species=s.get("SPECIES"), status="ok",
        pr_before=len(pr_b), pr_after=len(pr_a), hd_before=len(hd_b), hd_after=len(hd_a),
        pr_conf=";".join(x["confidence"] for x in pr_a), pr_class=";".join(x["locus_class"] for x in pr_a),
        pr_at_flagged=at_flag, flagged_loci=len(fl), pr_on_hd_contig=shared,
        receptors=nrec, precursors=npre,
        pr_withheld_after=sum(x["family"] == "Basidiomycota:PR" for x in (ra.get("suppressed_loci") or [])),
        wall_before=rb.get("_wall") if rb else None, wall_after=ra.get("_wall")))

keys = list(rows[0].keys()) if rows else []
allk = sorted({k for r in rows for k in r}, key=lambda k: (keys + [k]).index(k) if k in keys else 99)
with open(f"{HERE}/per_genome.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=allk, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

ok = [r for r in rows if r["status"] == "ok"]
for order in ("Polyporales", "Russulales"):
    rs = [r for r in ok if r["order"] == order]
    if not rs:
        continue
    called_b = sum(r["pr_before"] > 0 for r in rs)
    called_a = sum(r["pr_after"] > 0 for r in rs)
    print(f"{order}: {len(rs)} genomes | PR called before {called_b}, after {called_a} | "
          f"HD called before {sum(r['hd_before']>0 for r in rs)}, after {sum(r['hd_after']>0 for r in rs)}")
    pc = [r for r in rs if r["pr_after"] > 0]
    print(f"  PR calls {sum(r['pr_after'] for r in rs)}; at a flagged locus {sum(r['pr_at_flagged'] for r in rs)}; "
          f"genomes with flagged loci {sum(r['flagged_loci']>0 for r in rs)}; PR on an HD contig {sum(r['pr_on_hd_contig'] for r in rs)}")
    if pc:
        print(f"  receptors per called genome median {statistics.median(r['receptors'] for r in pc)}, "
              f"precursors median {statistics.median(r['precursors'] for r in pc)}")
        conf = {}
        for r in pc:
            for c in r["pr_conf"].split(";"):
                conf[c] = conf.get(c, 0) + 1
        print(f"  PR confidence {conf}")
    wb = [r["wall_before"] for r in rs if r["wall_before"]]
    wa = [r["wall_after"] for r in rs if r["wall_after"]]
    if wb and wa:
        print(f"  wall median before {statistics.median(wb):.0f} s, after {statistics.median(wa):.0f} s")
print("missing reports:", [r["asm"] for r in rows if r["status"] != "ok"])
