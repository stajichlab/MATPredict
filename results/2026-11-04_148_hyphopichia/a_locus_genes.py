import yaml,sys
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-04_148_hyphopichia/"
for g in ("GCA_030444895.1_ASM3044489v1","GCA_030556945.1_ASM3055694v1"):
    print("==",g)
    for run in ("base","cand"):
        r=yaml.safe_load(open(f"{O}{run}/runs/{g}/detection_report.yaml"))
        for l in r["detected"]:
            if l["idiomorph"] in ("A",):
                print(run,l["confidence"],l["contig"],l["start"],l["end"])
                for e in l["gene_evidence"]: print("    ",e["gene"],e["role"],round(e["identity"] or 0),e["coverage"] and round(e["coverage"]),round(e["bitscore"] or 0),e["status"],e["reference_record"][:22],e["method"])
        print("   suppressed:",[(s["family"],s["idiomorph"],s["genes_found"],s["withheld_reason"]) for s in r.get("suppressed_loci",[]) if s["contig"]==l["contig"]][:4])
