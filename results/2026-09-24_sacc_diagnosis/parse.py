import yaml,json,os,sys
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-23_saccharomyces_bar/runs/"
out=[]
for a in sorted(os.listdir(R)):
    p=R+a+"/detection_report.yaml"
    if not os.path.exists(p): out.append({"asm":a,"missing":True}); continue
    d=yaml.load(open(p),Loader=yaml.CSafeLoader)
    loci=[]
    for L in d.get("detected") or []:
        loci.append({k:L.get(k) for k in ["contig","start","end","confidence","idiomorph","idiomorph_candidates","locus_class","polished_genes","genes_found","idiomorph_margin","idiomorph_resolutions","fragmented","reference_records","detection_pass"]} | {"edge":[s.get("contig_edge_distance") for s in L["segments"]],
          "ev":[(e["gene"],e["start"],e["end"],e["strand"],e["identity"],e["status"],e["method"],e["reference_record"]) for e in L["gene_evidence"]]})
    ev=[];res=[];ps=[]
    dp=R+a+"/evidence_diagnostics.jsonl"
    if os.path.exists(dp):
        for l in open(dp):
            x=json.loads(l)
            if x["kind"]=="evidence": ev.append({k:x[k] for k in ["contig","cluster_start","cluster_end","gene_count","hit_count","roles","best_identity","admitted"]})
            elif x["kind"]=="idiomorph_resolution": res.append({k:x[k] for k in ["contig","winner","loser","winner_identity","loser_identity"]})
            else: ps.append({k:x[k] for k in ["contig","cluster_start","evidenced_idiomorph","rescues_skipped"]})
    out.append({"asm":a,"routing":d.get("routing_mode"),"suppressed":d.get("suppressed_unpolished"),"loci":loci,
      "not_detected":d.get("not_detected"),"evidence":ev,"resolutions":res,"polish_scope":ps})
json.dump(out,open("parsed.json","w"))
print(len(out))
