import json,collections,csv
D=json.load(open("parsed.json"))
C=collections.Counter
print("loci dist",C(len(g["loci"]) for g in D))
print("suppressed dist",sorted(C(g["suppressed"] for g in D).items())[:30])
# not_detected reasons
nd=C()
for g in D:
    for n in g["not_detected"] or []:
        r=n["reason"]; nd[r.split(",")[0][:90]]+=1
print(nd.most_common(10))
# genes_found signature per locus
sig=C(); idio=C()
for g in D:
    for L in g["loci"]:
        sig[(L["idiomorph"],tuple(sorted(L["genes_found"])))]+=1
for k,v in sig.most_common(20): print(v,k)
print(C(L["polished_genes"] for g in D for L in g["loci"]))
print(C(L["detection_pass"] for g in D for L in g["loci"]))
print(C(L["confidence"] for g in D for L in g["loci"]))
