import yaml,glob
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"
for g in ("GCA_051599425.1_Sai_col_YB2330_v1","GCF_001661265.1_Saico1"):
    f=glob.glob(f"{R}/2026-10-26_v070_ascomycota/wave_*/small_*/runs/{g}/detection_report.yaml")[0]
    r=yaml.safe_load(open(f)); print("==",g)
    L=r.get("suppressed_loci") or []
    print("suppressed loci:",len(L),"| single-gene:",sum(1 for l in L if len(l.get("genes_found") or [])==1))
    for l in L:
        if len(l.get("genes_found") or [])>=2: print("  multi:",{k:l.get(k) for k in ("contig","start","end","idiomorph","genes_found","polished_genes","withheld_reason")})
    print("receptor/other keys:",{k:r.get(k) for k in ("assembly_gap_at_locus","zygosity","two_idiomorphs")})
