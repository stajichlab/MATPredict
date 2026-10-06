import json,collections,csv
D=json.load(open("parsed.json"))
C=collections.Counter
res=C(); n=0
for g in D:
    L=g["loci"]
    if len(L)==2 and L[0]["contig"]==L[1]["contig"]:
        L=sorted(L,key=lambda l:l["start"]); d=L[1]["start"]-L[0]["end"]
        if 240e3<d<340e3:
            n+=1
            mid=[e for e in g["evidence"] if e["contig"]==L[0]["contig"] and L[0]["end"]+100e3<e["cluster_start"]<L[1]["start"]-40e3]
            best=max([e["best_identity"] for e in mid] or [0])
            res[int(best//10*10)]+=1
            if n<=4:
                print(g["asm"],L[0]["contig"],[(l["start"],l["idiomorph"]) for l in L]); 
                for e in mid: print("   ",e)
print(n,res)
