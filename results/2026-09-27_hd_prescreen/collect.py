"""Step 1: every redHD candidate cluster per genome, labelled positive (overlaps
the called HD locus) or negative, with its polish status from the run."""
import json, glob, os, yaml, csv
RUN = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_rhodotorula_curation/run_02434ee/runs"
FAM = "Basidiomycota:redHD"
rows = []
for d in sorted(glob.glob(f"{RUN}/*/")):
    g = os.path.basename(d.rstrip("/"))
    rep = yaml.safe_load(open(d + "detection_report.yaml"))
    called = [(x["contig"], x["start"], x["end"]) for x in rep.get("detected") or [] if x["family"] == FAM]
    supp = [(x["contig"], x["start"], x["end"]) for x in rep.get("suppressed_loci") or [] if x["family"] == FAM]
    wall = open(d + "wall_seconds").read().strip() if os.path.exists(d + "wall_seconds") else ""
    polished = {}
    for l in open(d + "evidence_diagnostics.jsonl"):
        r = json.loads(l)
        if r.get("family") != FAM:
            continue
        if r["kind"] == "polish_scope":
            polished[(r["contig"], r["cluster_start"], r["cluster_end"])] = r.get("polish_pairs_attempted")
    for l in open(d + "evidence_diagnostics.jsonl"):
        r = json.loads(l)
        if r.get("family") != FAM or r["kind"] != "evidence":
            continue
        c, s, e = r["contig"], r["cluster_start"], r["cluster_end"]
        ov = lambda L: any(cc == c and s <= ee and e >= ss for cc, ss, ee in L)
        rows.append(dict(genome=g, contig=c, start=s, end=e, gene_count=r["gene_count"],
                         hit_count=r["hit_count"], best_identity=r["best_identity"],
                         admitted=r["admitted"], polish_capped=r["polish_capped"],
                         polished=(c, s, e) in polished, label="pos" if ov(called) else "neg",
                         withheld=ov(supp), wall=wall))
with open("clusters.tsv", "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
import collections
c = collections.Counter((r["label"], r["admitted"], r["polish_capped"]) for r in rows)
print(len(rows), "clusters in", len({r['genome'] for r in rows}), "genomes"); print(c)
print("genomes with a pos:", len({r['genome'] for r in rows if r['label']=='pos'}))
