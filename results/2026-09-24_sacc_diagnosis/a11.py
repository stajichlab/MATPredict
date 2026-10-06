import json,collections,re,csv
S=json.load(open("sites.json")); D={g["asm"]:g for g in json.load(open("parsed.json"))}
M={r["ASMID"]:r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
C=collections.Counter(); ex=collections.defaultdict(list)
for a,s in S.items():
    for c,(st,l,rest) in s.items():
        if l!="not_called": continue
        m=re.match(r"gap=(\d+);N=(\d+);contig=(\S+);pos=(\d+)",rest); gap,N,ctg,pos=int(m[1]),int(m[2]),m[3],int(m[4])
        ev=[e for e in D[a]["evidence"] if e["contig"]==ctg and e["cluster_start"]<=pos+gap+500 and e["cluster_end"]>=pos-500]
        best=max(ev,key=lambda e:e["best_identity"],default=None)
        if best is None or best["best_identity"]<90: k="no >=90% cluster"
        else: k=f"cluster gene_count={best['gene_count']} admitted={best['admitted']}"
        C[(c,k)]+=1
        if len(ex[k])<3: ex[k].append((a,c,rest,best))
for k,v in sorted(C.items()): print(k,v)
for k,v in ex.items():
    for x in v: print("EX",k,x)
