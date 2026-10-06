import json,collections,csv,re
D={g["asm"]:g for g in json.load(open("parsed.json"))}
M={r["ASMID"]:r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
C=collections.Counter
exp={"HML":4499,"MAT":4499,"HMR":3499}
tab=C(); ex=collections.defaultdict(list); sp=C()
pergenome={}
for l in open("flank_all.tsv"):
    f=l.rstrip("\n").split("\t"); a=f[0]; g=D[a]
    species=M.get(a,{}).get("SPECIES","?"); cer = species=="Saccharomyces cerevisiae"
    sites={}
    for x in f[1:]:
        c,rest=x.split(":",1)
        m=re.match(r"gap=(\d+);N=(\d+);contig=(\S+);pos=(\d+)",rest)
        if m:
            gap,N,ctg,pos=int(m[1]),int(m[2]),m[3],int(m[4])
            if N>0: st="N_gap"
            elif gap>=exp[c]-600 and gap<=exp[c]+800: st="intact"
            elif gap<exp[c]-600: st="shortened"
            else: st="expanded"
            called=[L for L in g["loci"] if L["contig"]==ctg and L["start"]<=pos+gap+2000 and L["end"]>=pos-2000]
            lab=",".join(sorted(L["idiomorph"] for L in called)) or "not_called"
        else:
            st=rest.split("(")[0]; lab="n/a"; gap=N=None
        sites[c]=(st,lab,rest)
        grp="cerevisiae" if cer else "other"
        tab[(grp,c,st,"called" if lab not in("not_called","n/a") else lab)]+=1
        if len(ex[(c,st,lab=="not_called")])<3 and cer: ex[(c,st,lab=="not_called")].append((a,rest,lab))
    pergenome[a]=sites
json.dump(pergenome,open("sites.json","w"))
for k,v in sorted(tab.items()): print(k,v)
for k,v in ex.items():
    if k[2]:
        for x in v: print("EX",k,x)
