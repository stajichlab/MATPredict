import yaml,glob
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"
f=glob.glob(f"{R}/2026-10-26_v070_ascomycota/wave_*/small_*/runs/GCA_001929475.1_NeoIirr1.0/detection_report.yaml")[0]
r=yaml.safe_load(open(f))
print("families_attempted:",r.get("families_attempted"),"| not_detected:",str(r.get("not_detected"))[:300])
print("params:",{k:r["run"]["parameters"].get(k) for k in ("min_hits","min_identity","require_core_role","polish_strong_identity","polish_core_pair","max_polished_clusters_per_family")})
for l in r.get("suppressed_loci") or []:
    if l.get("contig") in ("LXFE01003167.1","LXFE01000411.1","LXFE01003595.1"):
        print(yaml.safe_dump(l,default_flow_style=True,width=250)[:900])
