import json,collections,csv,statistics as st
D=json.load(open("parsed.json"))
S={r["asm"]:r for r in csv.DictReader(open("asmstats.tsv"),delimiter="\t")}
C=collections.Counter
def ov(e,L): return e["contig"]==L["contig"] and e["cluster_start"]<=L["end"]+2000 and e["cluster_end"]>=L["start"]-2000
unrep=C(); ex=collections.defaultdict(list); gcount=C()
tot_unrep=0; genomes_with=0
for g in D:
    u=[e for e in g["evidence"] if e["best_identity"]>=90 and not any(ov(e,L) for L in g["loci"])]
    # dedupe clusters by coords
    u={(e["contig"],e["cluster_start"],e["cluster_end"]):e for e in u}.values()
    if u: genomes_with+=1
    for e in u:
        tot_unrep+=1; gcount[(e["gene_count"],e["admitted"])]+=1
        if len(ex[e["gene_count"]])<5: ex[e["gene_count"]].append((g["asm"],e))
print("genomes with >=1 unreported >=90% identity cluster:",genomes_with,"clusters:",tot_unrep)
print(gcount)
for k,v in ex.items():
    for x in v: print(k,x)
