import yaml,glob
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"
for g,c in (("GCA_051599425.1_Sai_col_YB2330_v1","AP039747.1"),("GCF_001661265.1_Saico1","NW_017567017.1")):
    f=glob.glob(f"{R}/2026-10-26_v070_ascomycota/wave_*/small_*/runs/{g}/detection_report.yaml")[0]
    r=yaml.safe_load(open(f))
    for l in r.get("suppressed_loci") or []:
        if l["contig"]==c and len(l.get("genes_found") or [])>=2:
            print("==",g); print(yaml.safe_dump({k:v for k,v in l.items()},default_flow_style=False,width=200)[:2500])
