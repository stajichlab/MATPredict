import yaml,glob,csv
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"
G=[r for r in csv.DictReader(open(R+"/2026-10-28_taphrinomycotina_mat_candidates/genomes.tsv"),delimiter="\t") if r["species"].startswith("Saitoella")]
for g in G:
    print("==",g["genome"],g["species"],g["status"])
    fs=glob.glob(f"{R}/2026-10-26_v070_ascomycota/wave_*/small_*/runs/{g['genome']}/detection_report.yaml")
    if not fs: print("no report"); continue
    r=yaml.safe_load(open(fs[0]))
    print("contigs",r["run"]["genome"].get("contigs"),"N50",r["run"]["genome"].get("n50"),"| not_detected:",str(r.get("not_detected"))[:350],"| suppressed:",r.get("suppressed_unpolished"))
    for l in (r.get("suppressed_loci") or [])[:6]:
        print("  ",{k:l.get(k) for k in ("contig","start","end","idiomorph","genes_found","polished_genes","withheld_reason")},"core_span",l.get("core_span"))
