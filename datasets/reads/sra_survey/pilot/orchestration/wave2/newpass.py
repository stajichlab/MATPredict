import csv,os,sys
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave1"
def rd(p): return list(csv.DictReader(open(p,newline="").read().replace("\r","").splitlines(),delimiter="\t"))
ok=("PASS","PROVISIONAL","PROVISIONAL_HOMOTHALLIC")
old={r["slug"]:r["decision"] for r in rd(W+"/gates.tsv")}; fin={r["slug"]:r["decision"] for r in rd(W+"/gates_final.tsv")}
new=[s for s,d in fin.items() if d in ok and old.get(s) not in ok]
print("passing in final:",sum(1 for d in fin.values() if d in ok),"| passed at test stage:",sum(1 for d in old.values() if d in ok),"| newly passing:",len(new))
done={f.split(".")[0] for f in os.listdir(W+"/out") if f.endswith(".run.json")}
tot=0;gb=0
for s in new:
    rows=rd(W+"/plan/"+s+"/runs_full.tsv")[:50]; todo=[r for r in rows if r["run_accession"] not in done]; tot+=len(todo); gb+=sum(int(r["base_count"] or 0) for r in todo)/1e9
    print(s,fin[s],"runs_full",len(rows),"to do",len(todo))
print("runs to do",tot,"Gbases",round(gb))
# also: passing species whose full stage is short of their cap
short=[]
for s,d in fin.items():
    if d in ok:
        rows=rd(W+"/plan/"+s+"/runs_full.tsv")[:50]; n=sum(1 for r in rows if r["run_accession"] not in done)
        if n and s not in new: short.append((s,n))
print("passing species with runs not yet run (not newly passing):",short)
