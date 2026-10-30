import csv,json,os,subprocess
P="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave2"
rd=lambda p:list(csv.DictReader(open(p,newline="").read().replace("\r","").splitlines(),delimiter="\t"))
rows=rd(f"{P}/plan/botrytis_cinerea/runs_full.tsv")
out=[]
for r in rows:
    run=r["run_accession"]; t=f"{P}/out/{run}.botrytis.tsv"; x=rd(t)[0] if os.path.exists(t) else {}
    out.append(dict(run=run,group="botrytis_cinerea",strain=r["strain"],study=r["study_accession"],call=x.get("call",""),shared_depth=x.get("shared_depth",""),truth=r["truth_detect"],sig=os.path.exists(f"{P}/out/{run}.sig.zip")))
os.makedirs(P+"/botrytis",exist_ok=True)
w=csv.DictWriter(open(P+"/botrytis/runs_calls.tsv","w"),list(out[0]),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(out)
print("with sketch:",sum(1 for o in out if o["sig"]),"of",len(out))
for o in out:
    if o["call"] not in ("MAT1-1","MAT1-2"): print(o["run"],o["strain"],o["study"],o["call"],o["shared_depth"],o["truth"])
