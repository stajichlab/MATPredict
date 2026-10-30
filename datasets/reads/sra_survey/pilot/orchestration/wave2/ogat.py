import yaml,glob,sys
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"
g=sys.argv[1]
f=glob.glob(f"{R}/2026-10-26_v070_ascomycota/wave_*/small_*/runs/{g}/detection_report.yaml")[0]
r=yaml.safe_load(open(f))
print("genome",g,"contigs",r["run"]["genome"].get("contigs"),"N50",r["run"]["genome"].get("n50"),"routing",r.get("routing_mode"),"| families_attempted",r.get("families_attempted"))
for l in r.get("detected") or []:
    print("\nDETECTED",{k:l.get(k) for k in ("family","contig","start","end","idiomorph","confidence","locus_class","detection_pass","polished_genes","high_basis")})
    print("   reference_records:",l.get("reference_records"))
    for e in l.get("gene_evidence") or []:
        print("   gene",e.get("gene"),e.get("role"),"| status",e.get("status"),"| id",round(e.get("identity",0),1) if e.get("identity") else None,"| bits",e.get("bitscore"),"| coords",e.get("start"),e.get("end"),"| record",e.get("reference_record_id"))
print("\nSUPPRESSED (withheld):")
for l in r.get("suppressed_loci") or []:
    if l.get("family") in ("Ascomycota:MTL","Ascomycota:MATsc","Ascomycota:MATog","Ascomycota:MATkp") and len(l.get("genes_found") or [])>=1:
        print("  ",{k:l.get(k) for k in ("family","contig","start","end","idiomorph","genes_found","polished_genes","withheld_reason")})
