import csv, collections, sys
D="data"
rd=lambda p:list(csv.DictReader(open(p),delimiter="\t"))
work={r[0]:r[2] for r in csv.reader(open(D+"/work.tsv"),delimiter="\t")}
gen={r["genome"]:r for r in rd("/Users/jstajich/projects/MATPredict/results/2026-10-03_ascomycota_v060/genomes.tsv")}
def load(p):
    d=collections.defaultdict(dict)
    for r in rd(f"{D}/loci_{p}.tsv"):
        if r["family"]=="Ascomycota:MAT": d[r["genome"]][(r["contig"],r["start"])]=r
    return d
L,C=load("legacy"),load("core_first")
S={(r["genome"],r["contig"],r["start"]):r for r in rd("ref/scored.tsv")}
# control diffs
n=0
for g in sorted(set(L)|set(C)):
    if work[g]!="ctl": continue
    for k in set(L[g])|set(C[g]):
        a,b=L[g].get(k),C[g].get(k)
        if (a or {}).get("idiomorph")!=(b or {}).get("idiomorph"):
            n+=1
            if n<=6: print(g,gen[g]["species"],k,"\n  L",a and (a["idiomorph"],a["genes_found"],a["cands"],a["n_res"]),"\n  C",b and (b["idiomorph"],b["genes_found"],b["cands"],b["n_res"]))
# changed loci vs profile
c=collections.Counter()
for g in L:
    if work[g]!="fix": continue
    for k in set(L[g])&set(C[g]):
        a,b=L[g][k],C[g][k]
        if a["idiomorph"]!=b["idiomorph"]:
            s=S.get((g,)+k)
            if s: c[(a["idiomorph"],b["idiomorph"],"old=prof" if a["idiomorph"]==s["profile_call"] else ("new=prof" if b["idiomorph"]==s["profile_call"] else "neither:"+s["profile_call"]))]+=1
for t,v in sorted(c.items()): print(t,v)
