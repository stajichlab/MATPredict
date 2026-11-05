"""Evidence per locus for the 80 genomes (run on HPCC): one row per detected locus plus suppressed loci."""
import glob, yaml, csv, json
O = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-04_189_secondary_loci/run/out/runs/"
rows = []
for f in sorted(glob.glob(O + "*/detection_report.yaml")):
    g = f.split("/")[-2]; d = yaml.safe_load(open(f))
    for l in d["detected"]:
        ge = [e for e in l["gene_evidence"] if e["role"] == "core_MAT" and not e.get("lost_to")]
        ids = [e["identity"] for e in ge if e["identity"] is not None]
        rows.append(dict(genome=g, kind="called", family=l["family"], contig=l["contig"], start=l["start"], end=l["end"], idiomorph=l["idiomorph"], confidence=l["confidence"],
                         locus_class=l["locus_class"], genes=",".join(sorted({e["gene"] for e in l["gene_evidence"] if not e.get("lost_to")})), core_identity_max=round(max(ids), 1) if ids else "",
                         core_identity_min=round(min(ids), 1) if ids else "", polished=l.get("polished_genes"), flags=",".join(l.get("resolution_flags", [])),
                         edge=min([s.get("contig_edge_distance") for s in l["segments"] if s.get("contig_edge_distance") is not None] or [""]), ambiguous=json.dumps(l.get("cross_call") or "")[:80]))
    for s in d.get("suppressed_loci", []):
        rows.append(dict(genome=g, kind="withheld:" + str(s.get("withheld_reason")), family=s["family"], contig=s["contig"], start=s["start"], end=s["end"], idiomorph=s.get("idiomorph"), genes=",".join(s.get("genes_found") or [])))
keys = sorted({k for r in rows for k in r}, key=lambda k: ["genome","kind","family","contig","start","end","idiomorph","confidence","locus_class","genes","core_identity_max","core_identity_min","polished","flags","edge","ambiguous"].index(k) if k in ["genome","kind","family","contig","start","end","idiomorph","confidence","locus_class","genes","core_identity_max","core_identity_min","polished","flags","edge","ambiguous"] else 99)
w = csv.DictWriter(open(O + "../evidence.tsv", "w"), fieldnames=keys, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows); print(len(rows))
