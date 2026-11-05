import yaml
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-04_148_hyphopichia/"
for g in ("GCA_030557965.1_","GCA_030566005.1_","GCA_030571075.1_","GCA_030578995.1_","GCA_030580355.1_","GCF_019202705.1_"):
    import glob
    d=glob.glob(O+"cand/runs/"+g+"*")[0]; r=yaml.safe_load(open(d+"/detection_report.yaml"))
    print("==",d.split("/")[-1][:28])
    for l in r["detected"]:
        if l["idiomorph"]=="alpha" and any(e["gene"] in("MTLA1","MTLA2") for e in l["gene_evidence"]):
            print("  alpha locus",l["contig"],l["start"],l["end"],l["confidence"])
            for e in sorted(l["gene_evidence"],key=lambda e:e["start"]):
                if e["role"]=="core_MAT": print("     ",e["gene"],e["start"],e["end"],e["strand"],"id",round(e["identity"] or 0),"cov",e["coverage"] and round(e["coverage"]),"bits",round(e["bitscore"] or 0),e["status"],e["reference_record"][:20])
