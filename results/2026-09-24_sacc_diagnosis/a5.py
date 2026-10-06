import json,collections
D=json.load(open("parsed.json"))
C=collections.Counter
def ov(e,L): return e["contig"]==L["contig"] and e["cluster_start"]<=L["end"]+2000 and e["cluster_end"]>=L["start"]-2000
h=C(); n=0; ng=C(); sumsup=0
for g in D:
    sumsup+=g["suppressed"] or 0
    for e in g["evidence"]:
        if e["admitted"] and not any(ov(e,L) for L in g["loci"]):
            n+=1; h[int(e["best_identity"]//10*10)]+=1; ng[e["gene_count"]]+=1
print("admitted unreported clusters",n,sorted(h.items()),ng, "sum suppressed",sumsup)
# all evidence identity hist
print(sorted(C(int(e["best_identity"]//10*10) for g in D for e in g["evidence"]).items()))
