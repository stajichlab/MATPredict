import json,collections,csv
D=json.load(open("parsed.json"))
S={r["asm"]:r for r in csv.DictReader(open("asmstats.tsv"),delimiter="\t")}
meta={r["ASMID"]:r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
C=collections.Counter
# per genome: high-identity evidence clusters (best_identity>=90)
def hi(g,t=90): return [e for e in g["evidence"] if e["best_identity"]>=t]
for n in [0,1,2,3]:
    gs=[g for g in D if len(g["loci"])==n]
    hc=C(len(hi(g)) for g in gs)
    print(n,len(gs),"hi-id(>=90) clusters per genome:",sorted(hc.items()))
    print("   species:",C(meta.get(g["asm"],{}).get("SPECIES","?") for g in gs).most_common(6))
gs=[g for g in D if len(g["loci"])==0]
for g in gs[:8]:
    print(g["asm"],S[g["asm"]]["n_contigs"],S[g["asm"]]["N50"],g["not_detected"][0]["reason"][:40],g["not_detected"][0].get("genes_found"))
    for e in hi(g,60): print("   ",e)
