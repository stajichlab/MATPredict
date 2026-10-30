import yaml,os,glob
E="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-28_taphrinomycotina_mat_candidates/e2e"
for g in ("GCA_001929475.1_NeoIirr1.0","GCA_051599425.1_Sai_col_YB2330_v1","GCF_001661265.1_Saico1"):
    for db in ("base","taph"):
        f=f"{E}/{g}_{db}/detection_report.yaml"
        if not os.path.exists(f): print(g,db,"no report"); continue
        r=yaml.safe_load(open(f))
        det=[l for l in r.get("detected") or [] if l.get("family")=="Ascomycota:PM"]
        nd=[x for x in r.get("not_detected") or [] if x.get("family")=="Ascomycota:PM"]
        sup=[l for l in r.get("suppressed_loci") or [] if l.get("family")=="Ascomycota:PM"]
        multi=[(l["contig"],l["start"],l["end"],l.get("genes_found"),l.get("polished_genes")) for l in sup if len(l.get("genes_found") or [])>=2]
        one=[(l["contig"],l["start"],l["end"],l.get("genes_found"),l.get("polished_genes")) for l in sup if len(l.get("genes_found") or [])==1 and (l.get("polished_genes") or 0)>=1]
        print(f"== {g.split('_')[0]} {db}: PM detected={len(det)} | suppressed PM loci={len(sup)} (polished>=1 single-gene: {len(one)}) | multi-gene suppressed: {multi}")
        for l in det: print("   DETECTED",{k:l.get(k) for k in ("contig","start","end","idiomorph","confidence","locus_class","polished_genes")},[ (x.get("gene"),x.get("status"),round(x.get("identity",0),1) if x.get("identity") else None) for x in l.get("gene_evidence") or []])
        for x in nd: print("   not_detected:",x.get("reason"))
        for o in one: print("   single polished:",o)
