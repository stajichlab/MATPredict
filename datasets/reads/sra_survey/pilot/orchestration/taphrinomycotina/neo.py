import yaml,glob,sys
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"
f=glob.glob(f"{R}/2026-10-26_v070_ascomycota/wave_*/small_*/runs/GCA_001929475.1_NeoIirr1.0/detection_report.yaml")[0]
r=yaml.safe_load(open(f))
print({k:r.get(k) for k in ("routing_mode","status","suppressed_unpolished","suppressed_flank_carried","suppressed_mat_gene_gate","suppressed_paralog_class","uncalled_reason","not_called_reason","reason")})
print("top keys:",list(r)[:30])
for l in (r.get("suppressed_loci") or [])[:8]:
    print({k:l.get(k) for k in ("family","contig","start","end","idiomorph","polished_genes","genes_found","reason","suppress_reason","confidence")})
for k in r:
    if "reason" in k or "uncalled" in k: print(k, str(r[k])[:300])
print("genome:",r["run"]["genome"])
