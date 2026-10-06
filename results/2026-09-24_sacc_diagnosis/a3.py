import json,collections,csv,statistics as st
D=json.load(open("parsed.json"))
S={r["asm"]:r for r in csv.DictReader(open("asmstats.tsv"),delimiter="\t")}
C=collections.Counter
def med(x): return st.median(x) if x else None
print("nloci\tn\tmed_contigs\tmed_N50\tmed_total\tfrac_contigs>200\tmax_id_median")
for n in [0,1,2,3,"4+"]:
    gs=[g for g in D if (len(g["loci"])==n if n!="4+" else len(g["loci"])>=4)]
    nc=[int(S[g["asm"]]["n_contigs"]) for g in gs if S[g["asm"]]["n_contigs"]!="NA"]
    n50=[int(S[g["asm"]]["N50"]) for g in gs if S[g["asm"]]["N50"]!="NA"]
    tot=[int(S[g["asm"]]["total"]) for g in gs if S[g["asm"]]["total"]!="NA"]
    mx=[max([e["best_identity"] for e in g["evidence"]] or [0]) for g in gs]
    print(n,len(gs),med(nc),med(n50),med(tot),round(sum(x>200 for x in nc)/len(nc),2),med(mx),sep="\t")
# zero genomes max identity histogram
z=[g for g in D if not g["loci"]]
print(C(int(max([e["best_identity"] for e in g["evidence"]] or [0])//10*10) for g in z))
# bins by N50
bins=[(0,20000),(20000,50000),(50000,100000),(100000,300000),(300000,700000),(700000,10**8)]
print("N50bin\tn\tmean_loci\tfrac0\tfrac3+\tboth_a_alpha")
for lo,hi in bins:
    gs=[g for g in D if S[g["asm"]]["N50"]!="NA" and lo<=int(S[g["asm"]]["N50"])<hi]
    if not gs: continue
    both=sum(1 for g in gs if {"MATa","MATalpha"}<= {L["idiomorph"] for L in g["loci"]})
    print(f"{lo}-{hi}",len(gs),round(st.mean(len(g["loci"]) for g in gs),2),round(sum(not g["loci"] for g in gs)/len(gs),2),round(sum(len(g["loci"])>=3 for g in gs)/len(gs),2),round(both/len(gs),2),sep="\t")
