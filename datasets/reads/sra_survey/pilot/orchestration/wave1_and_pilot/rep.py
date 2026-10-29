import sys,yaml,glob
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"
for g in sys.argv[1:]:
    d=glob.glob(f"{R}/2026-10-26_v070_ascomycota/wave_*/small_*/runs/{g}/detection_report.yaml")
    r=yaml.safe_load(open(d[0])); print("==",g,"| top keys:",[k for k in r if k not in("run",)][:14])
    for k in ("loci","detected_loci","calls"):
        for l in (r.get(k) or []):
            print(" ",{x:l.get(x) for x in ("family","contig","start","end","idiomorph","confidence","tier","locus_class","margin","genes_found","polished_genes","verification") if x in l})
    print("  withheld:",[(l.get("idiomorph"),l.get("contig"),l.get("genes_found")) for l in (r.get("suppressed_loci") or [])][:3])
