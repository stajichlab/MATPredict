import csv,os,collections,glob,json
P="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing"
rd=lambda p:list(csv.DictReader(open(p,newline="").read().replace("\r","").splitlines(),delimiter="\t"))
W=P+"/waves/wave2"
for stage in ("runs_test.tsv","runs_full.tsv"):
    rows=rd(f"{W}/plan/botrytis_cinerea/{stage}"); print(stage,len(rows),"runs planned")
rows=rd(f"{W}/plan/botrytis_cinerea/runs_full.tsv")
res=[]
for r in rows:
    run=r["run_accession"]; j=f"{W}/out/{run}.run.json"; t=f"{W}/out/{run}.botrytis.tsv"
    st=json.load(open(j))["state"] if os.path.exists(j) else "not_run"
    x=rd(t)[0] if os.path.exists(t) else {}
    res.append(dict(run=run,study=r["study_accession"],strain=r["strain"],layout=r["library_layout"],truth=r["truth_detect"],state=st,call=x.get("call",""),shared=x.get("shared_depth",""),flags=x.get("flags","")))
done=[x for x in res if x["state"]=="done"]
print("states:",dict(collections.Counter(x["state"] for x in res)))
c=collections.Counter(x["call"] for x in done); print("calls:",dict(c))
typed=[x for x in done if x["call"] in ("MAT1-1","MAT1-2")]
print("typed single-allele:",len(typed),"of",len(done),"finished;",{k:sum(1 for x in typed if x["call"]==k) for k in ("MAT1-1","MAT1-2")})
tr=[x for x in typed if x["truth"] in ("MAT1-1","MAT1-2")]; print("with assembly call:",len(tr),"agree:",sum(1 for x in tr if x["call"]==x["truth"]))
print("studies:",len({x['study'] for x in done}),dict(collections.Counter(x["study"] for x in done).most_common(8)))
print("both/none/low rows:",[(x["run"],x["strain"],x["call"],x["shared"],x["flags"]) for x in done if x["call"] not in ("MAT1-1","MAT1-2")])
# earlier pilots
for d in glob.glob(P+"/*/summary/summary_by_group.tsv")+glob.glob(P+"/pilot*/*/summary_by_group.tsv"):
    for r in rd(d):
        if "botry" in r.get("group","").lower(): print("earlier",d.split("/")[-3],r)
