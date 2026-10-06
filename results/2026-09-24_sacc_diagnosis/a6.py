import json,collections,csv
D=json.load(open("parsed.json"))
S={r["asm"]:r for r in csv.DictReader(open("asmstats.tsv"),delimiter="\t")}
C=collections.Counter
spans=[L["end"]-L["start"] for g in D for L in g["loci"]]
print("span max",max(spans),"n>5kb",sum(s>5000 for s in spans),"n>3kb",sum(s>3000 for s in spans))
ab={"MATa":"a","MATalpha":"A","undetermined":"u"}
pat3=C(); pat2=C(); samec=C(); ex=collections.defaultdict(list)
for g in D:
    L=g["loci"]
    contigs=C(l["contig"] for l in L)
    if len(L)==3 and len(contigs)==1:
        L=sorted(L,key=lambda l:l["start"])
        g1=L[1]["start"]-L[0]["end"]; g2=L[2]["start"]-L[1]["end"]
        if g1<g2: L=L[::-1]; g1,g2=g2,g1
        p="".join(ab[l["idiomorph"]] for l in L)
        pat3[p]+=1
        if len(ex[p])<3: ex[p].append((g["asm"],[(l["contig"],l["start"],l["idiomorph"],l["genes_found"]) for l in L],g1,g2))
    samec[(len(L),len(contigs))]+=1
print("3 loci on one contig, oriented HML..MAT..HMR (larger gap first):",pat3)
print("(n loci, n contigs):",sorted(samec.items()))
for p,v in ex.items():
    for x in v: print(p,x)
