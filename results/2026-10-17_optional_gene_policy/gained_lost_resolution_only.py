import csv, collections
rd=lambda p:list(csv.DictReader(open(p),delimiter="\t"))
sets=collections.defaultdict(set)
for r in csv.reader(open("data/work2.tsv"),delimiter="\t"):
    for s in r[2].split(","): sets[s].add(r[0])
gen={r["genome"]:r for r in rd("/Users/jstajich/projects/MATPredict/results/2026-10-03_ascomycota_v060/genomes.tsv")}
def load(p):
    return {(r["genome"],r["contig"],r["start"]):r for r in rd(p) if r["family"]=="Ascomycota:MAT"}
L=load("data/loci2_legacy.tsv"); R=load("data/loci2_resolution_only.tsv"); C=load("data/loci_core_first.tsv")
gained=[k for k in R if k not in L]; lost=[k for k in L if k not in R]
# pairing
def near(k,pool):
    return any(p[0]==k[0] and p[1]==k[1] and abs(int(p[2])-int(k[2]))<30000 for p in pool)
print("gained",len(gained),"lost",len(lost),"lost with a gained locus nearby (coordinate shift):",sum(near(k,gained) for k in lost))
c=collections.Counter()
for k in gained:
    r=R[k]; c[(r["confidence"],r["locus_class"])]+=1
print("gained conf/class",c.most_common())
print("gained gene sets",collections.Counter(R[k]["genes_found"] for k in gained).most_common(6))
print("set4 gained by order", collections.Counter((gen[k[0]].get("order") or gen[k[0]].get("class_")) for k in gained if k[0] in sets["set4"]).most_common(8))
# vs core_first gains for same genomes
inC=[k for k in gained if k in C]; print("gained under res_only and (genome in earlier run) also gained under core_first:", sum(1 for k in gained if k[0] in {x[0] for x in C} and k in C), "of", sum(1 for k in gained if k[0] in {x[0] for x in C}))
cg=[k for k in C if k not in {} ]
print("--- lost")
print(collections.Counter((L[k]["idiomorph"],L[k]["confidence"],L[k]["locus_class"],L[k]["genes_found"]) for k in lost).most_common(6))
print(collections.Counter(gen[k[0]]["species"].split()[0] for k in lost).most_common(6))
for k in lost[:2]:
    print(k, "genome other loci in res_only:", [(x[1],x[2],r["idiomorph"],r["genes_found"]) for x,r in R.items() if x[0]==k[0]])
rg={x[0] for x in R}
print("lost in genomes still having a locus:", sum(1 for k in lost if k[0] in rg), "of", len(lost))
print("gained loci in genomes that had no locus in legacy:", sum(1 for k in gained if not any(x[0]==k[0] for x in L)), "genomes:", len({k[0] for k in gained if not any(x[0]==k[0] for x in L)}))
