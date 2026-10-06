import json,collections,csv
D=json.load(open("parsed.json"))
S={r["asm"]:r for r in csv.DictReader(open("asmstats.tsv"),delimiter="\t")}
C=collections.Counter
ab={"MATa":"a","MATalpha":"A","undetermined":"u"}
def gapclass(d):
    if 60e3<d<140e3: return "MAT-HMR(~94kb)"
    if 150e3<d<240e3: return "HML-MAT(~187kb)"
    if 240e3<d<340e3: return "HML-HMR(~280kb)"
    return f"other({int(d/1000)}kb)"
two=C(); ex=collections.defaultdict(list)
for g in D:
    L=g["loci"]
    if len(L)==2 and L[0]["contig"]==L[1]["contig"]:
        L=sorted(L,key=lambda l:l["start"]); d=L[1]["start"]-L[0]["end"]
        k=gapclass(d); two[k]+=1
        if len(ex[k])<3: ex[k].append((g["asm"],S[g["asm"]]["n_contigs"],[(l["contig"],l["start"],l["end"],l["idiomorph"],l["edge"]) for l in L]))
print("2 loci same contig:",two)
for k,v in ex.items():
    for x in v: print(k,x)
# edge distance of all reported loci
ed=[min(e for e in l["edge"] if e is not None) for g in D for l in g["loci"] if any(e is not None for e in l["edge"])]
print("reported loci with contig_edge_distance <1kb:",sum(e<1000 for e in ed),"<5kb",sum(e<5000 for e in ed),"of",len(ed))
