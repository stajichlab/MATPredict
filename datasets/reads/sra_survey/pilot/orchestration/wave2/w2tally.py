import os,csv,collections
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave2/plan"
cat=collections.defaultdict(list)
for s in sorted(os.listdir(W)):
    st=open(f"{W}/{s}/status").read().strip() if os.path.exists(f"{W}/{s}/status") else "planned"
    if st=="planned": cat["planned"].append(s); continue
    lt=[f for f in os.listdir(f"{W}/{s}/panel") if f.endswith("_lineages.tsv")] if os.path.isdir(f"{W}/{s}/panel") else []
    if not lt: cat[st+" | no loci at all (no MAT1-1/MAT1-2 loci found)"].append(s); continue
    rows=list(csv.DictReader(open(f"{W}/{s}/panel/{lt[0]}"),delimiter="\t"))
    a=sum(int(r["MAT1-1"]) for r in rows); b=sum(int(r["MAT1-2"]) for r in rows)
    k=f"{st} | "+("no MAT1-1/MAT1-2 loci (other idiomorph names or none)" if a==0 and b==0 else "only MAT1-1 in assemblies" if b==0 and a else "only MAT1-2 in assemblies" if a==0 and b else "both idiomorphs")
    cat[k].append(s)
for k,v in cat.items(): print(len(v),k,"->",", ".join(v[:6]))
