import json,collections,re,csv
S=json.load(open("sites.json"))
M={r["ASMID"]:r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
C=collections.Counter
lab=C(); nb=C()
for a,s in S.items():
    if M.get(a,{}).get("SPECIES")!="Saccharomyces cerevisiae": continue
    for c,(st,l,rest) in s.items():
        if st=="intact": lab[(c,l)]+=1
        if st=="N_gap":
            m=re.match(r"gap=(\d+);N=(\d+)",rest); gap,N=int(m[1]),int(m[2]); nonN=gap-N
            b="nonN<1000" if nonN<1000 else ("1000-2500" if nonN<2500 else ">=2500")
            b2="N<=10" if N<=10 else ("N 11-500" if N<=500 else "N>500")
            nb[(c,b,l=="not_called")]+=1
            nb[("Nsize",c,b2,l=="not_called")]+=1
print("labels at intact sites:"); [print(" ",k,v) for k,v in sorted(lab.items())]
print("N-gap sites, (site, non-N bp, not_called):"); [print(" ",k,v) for k,v in sorted(nb.items(),key=str)]
